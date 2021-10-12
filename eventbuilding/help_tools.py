#!/usr/bin/env python
from __future__ import print_function
import os
import numpy as np


class EcalNumbers:
    def __init__(self):
        self.n_chips = 16
        self.n_scas = 15
        self.n_channels = 64
        self.n_slabs = 15  # Can be overwritten in CLI.
        self.cob_slabs = {5, 8}  # TODO: Check the COB positions.

        # TODO What are the real pos_z and Tungsten values for 15 slabs?
        self.pos_z = np.array([0, 2, 4, 6, 8, 9, 12, 14, 16, 18, 20, 22, 24, 26, 28]) * 15 # mm gap.
        self.w_config = {  # abs thickness of Tungsten/W plates.
            1: np.array([0, 2.1, 2.1, 4.2, 4.2, 0, 4.2, 2.1, 2.1, 0, 0, 0, 0, 0, 0]),
        }

        self.bcid_val_event = 50
        self.bcid_merge_delta = 30
        self.bcid_too_many_hits = 8000

        self.pedestal_min_average = 200
        self.pedestal_min_scas = 6
        self.pedestal_min_value = 10
        self.mip_cutoff = 0.5
        self.validate_ecal_numbers(self)


    @classmethod
    def validate_ecal_numbers(cls, n):
        assert type(n.n_chips) == int
        assert type(n.n_scas) == int
        assert type(n.n_channels) == int
        assert type(n.n_slabs) == int
        assert all((i_cob < n.n_slabs for i_cob in n.cob_slabs))

        assert type(n.pos_z) == np.ndarray
        assert len(n.pos_z) == n.n_slabs
        assert all(type(w_conf) == np.ndarray for w_conf in n.w_config.values())
        assert all(len(w_conf) == n.n_slabs for w_conf in n.w_config.values())

        assert type(n.bcid_val_event) == int
        assert n.bcid_merge_delta >= 0
        assert type(n.bcid_too_many_hits) == int

        assert type(n.pedestal_min_average) == int
        assert type(n.pedestal_min_scas) == int
        assert type(n.pedestal_min_value) == int
        assert type(n.mip_cutoff) == float and n.mip_cutoff <= 1


class EventBuildingException(Exception):
    pass


class EcalHit:
    def __init__(self, slab, chip, chan, sca, hg, lg, gain_hit_high, ecal_config):
        self.slab = slab
        self.chip = chip
        self.chan = chan
        self.sca = sca
        self.hg = hg
        self.lg = lg
        self._ecal_config = ecal_config

        self.isMasked = int(self._ecal_config.masked_map[self.slab][self.chip][self.chan])
        self.isCommissioned = 1
        self._gain_hit_high = gain_hit_high

        self._set_positions()
        self._pedestal_subtraction()
        self._mip_calibration()


    @property
    def isHit(self):
        return self._gain_hit_high


    def _set_positions(self):
        self.x0 = self._ecal_config.pos_x0[self.slab]
        self.z = self._ecal_config._N.pos_z[self.slab]
        slab_channel_map = self._ecal_config.get_channel_map(self.slab)
        (self.x,self.y) = slab_channel_map[(self.chip, self.chan)]
        # TODO: Is this really doing the right thing to the positioning? Check with new cosmic runs.
        # Invert the mapping for the 5 first slabs, to agree with the last 4.
        # if slab < 5: currently equivalent to:
        # if "_dif_" in self._ecal_config._N.slab_map[self.slab]:
        #    self.x = -self.x
        #    self.y = -self.y


    def _pedestal_subtraction(self):
        pedestals_per_sca = self._ecal_config.pedestal_map[self.slab][self.chip][self.chan]
        is_good_pedestal = pedestals_per_sca > self._ecal_config._N.pedestal_min_average
        if sum(is_good_pedestal) < self._ecal_config._N.pedestal_min_scas:
            self.isCommissioned = 0

        sca_has_valid_pedestal = pedestals_per_sca[self.sca] > self._ecal_config._N.pedestal_min_value
        if sca_has_valid_pedestal:
            self.hg -= pedestals_per_sca[self.sca]
        else:
            if (not self.isMasked) and self.isCommissioned:
                pedestal_average = np.mean(pedestals_per_sca[is_good_pedestal])
                self.hg -= pedestal_average
                self.isCommissioned = 0


    def _mip_calibration(self):
        mip_value = self._ecal_config.mip_map[self.slab][self.chip][self.chan]
        if mip_value > self._ecal_config._N.mip_cutoff:
            self.energy = self.hg / mip_value
        else:
            self.energy = 0
            self.isCommissioned = 0


class EcalConfig:

    def __init__(
        self,
        commissioning_tag,
        w_config=1,
        mapping_file="mapping/fev10_chip_channel_x_y_mapping.txt",
        mapping_file_cob="mapping/fev11_cob_chip_channel_x_y_mapping.txt",
        pedestals_file="pedestals/pedestal{tag}.txt",
        mip_calibration_file="mip_calib/MIP{tag}.txt",
        masked_file="masked/masked{tag}.txt",
        commissioning_folder=None,
        numbers=None,
        error_on_missing_config=True,
        verbose=False,
    ):
        self._commissioning_tag = commissioning_tag
        self._verbose = verbose
        self._error_on_missing_config = error_on_missing_config
        if commissioning_folder:
            self._commissioning_folder = commissioning_folder
        else:
            # Resolves to the root folder of this repo
            # Equivalent to ../ if called from within the eventbuilding folder.
            self._commissioning_folder = os.path.dirname(os.path.dirname((__file__)))
        if numbers:
            EcalNumbers.validate_ecal_numbers(numbers)  # Catch problems early on.
            self._N = numbers
        else:
            self._N = EcalNumbers()

        self.pos_x0 = self._build_w_config(w_config)
        self._channel_map = self._read_mapping(mapping_file)
        self._channel_map_cob = self._read_mapping(mapping_file_cob)
        self.pedestal_map = self._read_pedestals(pedestals_file)
        self.mip_map = self._read_mip_values(mip_calibration_file)
        self.masked_map = self._read_masked(masked_file)


    def get_channel_map(self, slab):
        if slab in self._N.cob_slabs:
            return self._channel_map
        else:
            return self._channel_map_cob


    def _build_w_config(self, config):
        if config in self._N.w_config.keys():
            abs_thick = self._N.w_config[config]
        elif config == 0:
            abs_thick = np.zeros_like(self._N.pos_z)
        else:
            raise EventBuildingException("Not a valid W config:", config)

        w_x0 = 1 / 3.5  # 0.56 #Xo per mm of W.
        pos_x0 = np.cumsum(abs_thick) * w_x0
        print("W config %i used:" %config)
        print("absolute thickness:", abs_thick, "\npos_x0:", pos_x0)
        return pos_x0


    def _get_lines(self, file_name):
        if "{tag}" in file_name:
            file_name = file_name.format(tag=self._commissioning_tag)
        if not os.path.isabs(file_name):
            file_name = os.path.join(self._commissioning_folder, file_name)
        if not os.path.exists(file_name):
            raise EventBuildingException("File does not exist: %s" %file_name)

        try:
            lines = open(file_name).readlines()
        except FileNotFoundError:
            txt = "%s does not exist" %file_name
            if self._error_on_missing_config:
                raise EventBuildingException(txt)
        return lines


    def _read_mapping(self, file_name):
        channel_map = dict()
        lines = self._get_lines(file_name)

        fields = lines[0].split()
        i_chip = fields.index("chip")
        i_channel = fields.index("channel")
        i_x = fields.index("x")
        i_y = fields.index("y")

        for line in lines[1:]:
            v = line.split()
            pos = (float(v[i_x]), float(v[i_y]))
            channel_map[(int(v[i_chip]), int(v[i_channel]))] = pos
        return channel_map


    def _read_pedestals(self, file_name):
        ped_map = np.zeros((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
            self._N.n_scas,
        ))
        print("Reading pedestals from %s." %file_name)
        lines = self._get_lines(file_name)

        assert lines[0].startswith("#pedestal results")
        assert lines[1].startswith("#") and not lines[2].startswith("#")
        fields = lines[1][1:].split()
        i_slab = fields.index("layer")
        i_chip = fields.index("chip")
        i_channel = fields.index("channel")
        i_ped0 = fields.index("ped0")


        for line in lines[2:]:
            v = line.split()
            ped_val = float(v[i_ped0])  # TODO: Is it really ok to just us the pedestal mean from the first SCA?
            # TODO: Proposal:
            # for i_sca in range(ped_map.shape[-1]):
            #     n_entries_per_sca = 3  # ped, eped, widthped
            #     ped_val = float(v[2 + i_sca * n_entries_per_sca])
            #     ped_map[int(v[i_slab])][int(v[i_chip])][int(v[i_channel])][i_sca] = ped_val
            ped_map[int(v[i_slab])][int(v[i_chip])][int(v[i_channel])] = ped_val
        if self._verbose:
            print("pedestal_map", ped_map)
        # TODO: Are we really ok with 0-entries in the pedestal map?
        return ped_map


    def _read_mip_values(self, file_name):
        mip_map = np.ones((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
        ))
        print("Reading MIP values from %s." %file_name)
        lines = self._get_lines(file_name)

        assert lines[0].startswith("#mip results")
        assert lines[1].startswith("#") and not lines[2].startswith("#")
        fields = lines[1][1:].split()
        i_slab = fields.index("layer")
        i_chip = fields.index("chip")
        i_channel = fields.index("channel")
        i_mpv = fields.index("mpv")

        for line in lines[2:]:
            v = line.split()
            mip_val = float(v[i_mpv])
            mip_map[int(v[i_slab])][int(v[i_chip])][int(v[i_channel])] = mip_val
        if self._verbose:
            print("mip_map", mip_map)
        return mip_map


    def _read_masked(self, file_name):
        masked_map = np.zeros((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
        ))
        print("Reading masked channels from %s." %file_name)
        lines = self._get_lines(file_name)

        start_tag = "#masked_chns_list "
        assert lines[0].startswith(start_tag)
        assert not lines[1].startswith("#")
        fields = lines[0][len(start_tag):].split()
        assert fields[0] == "layer"
        assert fields[1] == "chip"
        assert fields[2] == "chns"

        for line in lines[1:]:
            v = line.split()
            assert len(v) == self._N.n_channels + 2
            slab = int(v[0])
            chip = int(v[1])
            for channel, mask_val in enumerate(v[2:]):
                masked_map[slab][chip][channel] = mask_val
        if self._verbose:
            print("masked_map", masked_map)
        return masked_map


if __name__ == "__main__":
    EcalConfig(verbose=True)
