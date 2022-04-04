import numpy as np


class BCIDHandler:
    def __init__(self, bcid_config, geometry_config, min_slabs_hit=1, merge_within_chip=False):
        self._skip_noisy_acquisition_start = int(bcid_config["skip_noisy_acquisition_start"])
        self._bad_value = int(bcid_config["bad_value"])
        self._drop_values = list(map(int, bcid_config["drop_values"].split(",")))
        self._drop_retrigger_delta = int(bcid_config["drop_retrigger_delta"])
        self._merge_delta = int(bcid_config["merge_delta"])
        self._overflow = int(bcid_config["overflow"])

        self._n_chips = int(geometry_config["n_chips"])
        self._n_scas = int(geometry_config["n_scas"])

        self._min_slabs_hit = min_slabs_hit
        self._merge_within_chip = merge_within_chip


    def _calculate_overflow_correction(self, arr):
        """"Within a chips memory, the BCID must increase. A drop indicates overflow."""
        # assert np.max(arr) < self._overflow
        n_memory_cycles = np.cumsum(arr[:,1:] - arr[:,:-1] < 0, axis=1)
        is_filled = arr[:,1:] != self._bad_value
        overflown_bcids = np.copy(arr)
        overflown_bcids[:,1:] = arr[:,1:] + self._overflow * n_memory_cycles * is_filled
        bcid_to_overflow_adjusted = {}
        unique_overflow = np.unique(overflown_bcids)
        unique_overflow = unique_overflow[unique_overflow >= self._overflow]
        for of_bcid in np.unique(overflown_bcids):
            bcid_to_overflow_adjusted[of_bcid % self._overflow] = of_bcid
        return overflown_bcids, bcid_to_overflow_adjusted


    def _is_retrigger(self, arr):
        """Within each chip's memory, look for consecutive BCIDs."""
        delta = arr[:,1:] - arr[:,:-1]
        return np.logical_and(delta > 0, delta <= self._drop_retrigger_delta)


    def _is_empty_sca(self, nhits_branch):
        """Drop the SKIROC2 empty SCAs."""
        nhits = np.array(list(nhits_branch)).reshape(-1, self._n_scas)
        return nhits == 0


    def _find_bcid_filling_last_sca(self, arr):
        bcid_in_last_sca = arr[:,-1]
        bcid_in_last_sca = bcid_in_last_sca[bcid_in_last_sca != self._bad_value]
        if len(bcid_in_last_sca):
            bcid_in_last_sca = bcid_in_last_sca.min()
        else:
            bcid_in_last_sca = self._bad_value
        return bcid_in_last_sca


    def _get_bcids(self, spill_entry):
        raw_bcids = np.array(list(spill_entry.bcid)).reshape(-1, self._n_scas)
        raw_bcids[self._is_empty_sca(spill_entry.nhits)] = self._bad_value
        chip_ids = np.nonzero(
            np.any(raw_bcids != self._bad_value, axis=1)
        )[0]
        bcids = np.copy(raw_bcids[chip_ids])
        if not self._merge_within_chip:
            bcids[:,1:][self._is_retrigger(bcids)] = self._bad_value
        bcids[bcids < self._skip_noisy_acquisition_start] = self._bad_value
        for drop_bcid in self._drop_values:
            bcids[bcids == drop_bcid] = self._bad_value
        return chip_ids, bcids


    def _get_bcid_start_stop(self, arr):
        all_bcids = np.unique(arr)
        all_bcids = all_bcids[all_bcids != self._bad_value]
        if len(all_bcids) == 0:
            return dict()
        previous_bcid = all_bcids[0]
        bcid_start_stop = {previous_bcid: previous_bcid}
        for current_bcid in all_bcids[1:]:
            if current_bcid - bcid_start_stop[previous_bcid] < self._merge_delta:
                bcid_start_stop[previous_bcid] = current_bcid
            else:
                bcid_start_stop[current_bcid] = current_bcid
                previous_bcid = current_bcid
        return bcid_start_stop


    def _get_array_positions_per_bcid(self, chip_ids, bcids):
        bcid_start_stop = self._get_bcid_start_stop(bcids)
        _pos_per_bcid = {k: {} for k in bcid_start_stop}
        bcid_before_merge = {k: {} for k in bcid_start_stop}
        for bcid_start in bcid_start_stop:
            i_chip_local, i_sca = np.nonzero(np.logical_and(
                bcids >= bcid_start, bcids <= bcid_start_stop[bcid_start]
            ))
            i_chip = chip_ids[i_chip_local]
            for i in range(len(i_chip_local)):
                i_chip = chip_ids[i_chip_local[i]]
                if i_chip in _pos_per_bcid[bcid_start]:
                    if self._merge_within_chip:
                        # Retrigger handling is postponed to build_events.py
                        # when merging within a chip.
                        _pos_per_bcid[bcid_start][i_chip].append(i_sca[i])
                    else:
                        # If we do not merge within a chip, then we take the first
                        # BCID that is written to the BCID within the event's BCID
                        # window and that is not an empty event (avoids retriggers).
                        continue
                else:
                    _pos_per_bcid[bcid_start][i_chip] = [i_sca[i]]
                bcid_before_merge[bcid_start][i_chip] = bcids[i_chip_local[i], i_sca[i]]

        pos_per_bcid = {}
        for bcid_start, chip_to_scas in _pos_per_bcid.items():
            chips = np.array(list(chip_to_scas.keys()))
            slabs_hit = chips // self._n_chips
            if len(np.unique(slabs_hit)) >= self._min_slabs_hit:
                pos_per_bcid[bcid_start] = chip_to_scas
        return pos_per_bcid, bcid_before_merge, bcid_start_stop


    def load_spill(self, spill_entry):
        """False if the current spill is empty."""
        chip_ids, bcids = self._get_bcids(spill_entry)
        overflown_bcids, self.bcid_to_overflow_adjusted = \
            self._calculate_overflow_correction(bcids)
        self.bcid_first_sca_full = self._find_bcid_filling_last_sca(overflown_bcids)
        self.pos_per_bcid, self._bcid_before_merge, bcid_start_stop = \
            self._get_array_positions_per_bcid(chip_ids, bcids)
        self.spill_bcids_overflown = []
        for start_bcid in self.pos_per_bcid.keys():
            n_memory_cycles = max([
                self.bcid_to_overflow_adjusted.get(b, b)
                for b in range(start_bcid, bcid_start_stop[start_bcid] + 1)
            ]) // self._overflow
            _overflown = start_bcid + n_memory_cycles * self._overflow
            self.spill_bcids_overflown.append(_overflown)
        self.spill_bcids_overflown = sorted(self.spill_bcids_overflown)
        return len(self.pos_per_bcid) > 0


    def yield_from_spill(self):
        for overflown_bcid in self.spill_bcids_overflown:
            bcid = overflown_bcid % self._overflow
            before_merge = self._bcid_before_merge[bcid]
            plus_overflow = (overflown_bcid // self._overflow) * self._overflow

            for chip, scas in self.pos_per_bcid[bcid].items():
                yield chip, scas, bcid + plus_overflow, before_merge[chip] + plus_overflow


    def previous_bcid(self, overflown_bcid):
        i_bcid = self.spill_bcids_overflown.index(overflown_bcid)
        if i_bcid == 0:
            return -1
        else:
            return self.spill_bcids_overflown[i_bcid - 1]

    def next_bcid(self, overflown_bcid):
        i_bcid = self.spill_bcids_overflown.index(overflown_bcid)
        if i_bcid == len(self.spill_bcids_overflown) - 1:
            return -1
        else:
            return self.spill_bcids_overflown[i_bcid + 1]
