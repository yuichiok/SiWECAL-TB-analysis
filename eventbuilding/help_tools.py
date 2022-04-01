#!/usr/bin/env python
from __future__ import print_function
import os
import sys

import numpy as np


asu_types = ["10", "11", "12" , "13", "COB"]


class EcalNumbers:
    def __init__(self, slabs=None, asu_version=None):
        self.n_chips = 16
        self.n_scas = 15
        self.n_channels = 64
        if slabs is None:
            self.slabs = []
        else:
            self.slabs = sorted(set(slabs))
        if asu_version is None or asu_version == "":
            asu_version = ["12" for s in self.slabs]
        elif type(asu_version) == str:
            asu_version = asu_version.split(",")
        self.asu_version = [asu_v.strip() for asu_v in asu_version]
        self.n_slabs = len(self.slabs)

        self.pedestal_min_average = 200
        self.pedestal_min_scas = 3
        self.pedestal_min_value = 10
        self.mip_cutoff = 0.5
        self.mip_malfunctioning_chip = 1000
        self.validate_ecal_numbers(self)


    @classmethod
    def validate_ecal_numbers(cls, n):
        assert type(n.n_chips) == int
        assert type(n.n_scas) == int
        assert type(n.n_channels) == int
        assert type(n.n_slabs) == int
        assert len(n.slabs) == n.n_slabs
        assert all((type(i_slab) == int for i_slab in n.slabs))
        assert len(n.asu_version) == len(n.slabs), n.asu_version + n.slabs
        assert all((asu_v in asu_types for asu_v in n.asu_version)), n.asu_version

        assert type(n.pedestal_min_average) == int
        assert type(n.pedestal_min_scas) == int
        assert type(n.pedestal_min_value) == int
        assert type(n.mip_cutoff) == float and n.mip_cutoff <= 1
        assert type(n.mip_malfunctioning_chip) in [int, float]  and n.mip_malfunctioning_chip != 0


class EventBuildingException(Exception):
    pass

dummy_config = dict(
    mapping_file="mapping/fev10_chip_channel_x_y_mapping.txt",
    mapping_file_cob="mapping/fev11_cob_chip_channel_x_y_mapping.txt",
    pedestals_file="pedestals/pedestal_PROTO15_dummy.txt",
    mip_calibration_file="mip_calib/MIP_PROTO15_dummy.txt",
    pedestals_lg_file="pedestals/pedestal_PROTO15_dummy_lg.txt",
    mip_calibration_lg_file="mip_calib/MIP_PROTO15_dummy_lg.txt",
    masked_file="masked/masked_PROTO15_dummy.txt",
)


def aligned_path(text, path):
    _text_length = 29
    assert _text_length >= len(text)
    return text + (_text_length - len(text)) * " " + str(path)


class EcalConfig:

    def __init__(
        self,
        calibration_files=dummy_config,
        numbers=None,
        zero_suppress=True,
        no_lg=False,
        error_on_missing_config=True,
        verbose=False,
    ):
        self._verbose = verbose
        self._error_on_missing_config = error_on_missing_config
        if numbers:
            EcalNumbers.validate_ecal_numbers(numbers)  # Catch problems early on.
            self._N = numbers
        else:
            self._N = EcalNumbers()

        self.x, self.y = self._get_x_y(
            calibration_files["mapping_file"],
            calibration_files["mapping_file_cob"],
        )
        self.z = self._get_z()

        self.masked_map = self._read_masked(calibration_files["masked_file"])
        self.pedestal_map = self._read_pedestals(calibration_files["pedestals_file"])
        self.mip_map = self._read_mip_values(calibration_files["mip_calibration_file"])

        self.is_commissioned_map = self._handle_uncommissioned_positions(
            self.pedestal_map, self.mip_map, self.masked_map
        )

        self.zero_suppress = zero_suppress
        self.no_lg = no_lg
        if self.no_lg:
            print("As requested with --no_lg, low gain will not be calibrated.")
        else:
            try:
                self.pedestal_lg_map = self._read_pedestals(calibration_files["pedestals_lg_file"])
                self.mip_lg_map = self._read_mip_values(calibration_files["mip_calibration_lg_file"])

                lg_is_commissioned = self._handle_uncommissioned_positions(
                    self.pedestal_lg_map, self.mip_lg_map, self.masked_map
                )
                self.is_commissioned_map = np.logical_and(
                    self.is_commissioned_map, lg_is_commissioned
                )
            except EventBuildingException as e:
                print("To run without the lowgain information, use --no_lg.")
                raise e


    def _get_x_y(self, mapping_file, mapping_file_cob):
        _x, _y = self._read_mapping_xy(mapping_file)
        _x_cob, _y_cob = self._read_mapping_xy(mapping_file_cob)

        xs, ys = [], []
        for slab in self._N.slabs:
            asu_version = self._N.asu_version[self._N.slabs.index(slab)]
            if asu_version == "COB":
                _xi = _x_cob
                _yi = _y_cob
            elif asu_version in ["10", "11", "12", "13"]:
                _xi = _x
                _yi = _y
            else:
                raise EventBuildingException(
                    "Asu version not recognized: " + str(asu_version)
                )
            if asu_version == "13":
                _xi = _xi + 60
            xs.append(_xi)
            ys.append(_yi)
        x = np.stack(xs)
        y = np.stack(ys)
        x, y = -x, -y
        return x, y


    def _get_z(self):
        """15 mm distance between slabs in the prototype."""
        def z_val_layer(slab):
            return np.full((self._N.n_chips, self._N.n_channels), 15 * slab)

        return np.stack([z_val_layer(slab) for slab in self._N.slabs])


    def _get_lines(self, file_name):
        try:
            lines = open(file_name).readlines()
        except FileNotFoundError:
            txt = "%s does not exist" %file_name
            if self._error_on_missing_config:
                raise EventBuildingException(txt)
        return lines


    def _read_mapping_xy(self, file_name):
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

        x = np.empty((self._N.n_chips, self._N.n_channels))
        y = np.empty((self._N.n_chips, self._N.n_channels))
        assert len(channel_map) == x.size
        for (chip, channel), (x_i, y_i) in channel_map.items():
            x[chip][channel] = x_i
            y[chip][channel] = y_i
        return x, y


    def _read_pedestals(self, file_name):
        ped_map = np.zeros((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
            self._N.n_scas,
        ))
        print(aligned_path("Reading pedestals from ", file_name))
        lines = self._get_lines(file_name)

        assert lines[0].startswith("#pedestal results")
        assert lines[1].startswith("#") and not lines[2].startswith("#")
        fields = lines[1][1:].split()
        i_slab = fields.index("layer")
        i_chip = fields.index("chip")
        i_channel = fields.index("channel")
        i_ped0 = fields.index("ped0")
        i_ped1 = fields.index("ped1")
        n_entries_per_sca = i_ped1 - i_ped0


        for line in lines[2:]:
            v = line.split()
            for i_sca in range(ped_map.shape[-1]):
                ped_val = float(v[i_ped0 + i_sca * n_entries_per_sca])
                err_ped_val = float(v[4 + i_sca * n_entries_per_sca])
                if err_ped_val < 0:
                    ped_val = 0
                idx_slab = self._N.slabs.index(int(v[i_slab]))
                ped_map[idx_slab][int(v[i_chip])][int(v[i_channel])][i_sca] = ped_val
        if self._verbose:
            print("pedestal_map", ped_map)
        return ped_map


    def _read_mip_values(self, file_name):
        mip_map = np.ones((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
        ))
        print(aligned_path("Reading MIP values from ", file_name))
        lines = self._get_lines(file_name)

        assert lines[0].startswith("#mip results")
        assert lines[1].startswith("#") and not lines[2].startswith("#")
        fields = lines[1][1:].split()
        i_slab = fields.index("layer")
        i_chip = fields.index("chip")
        i_channel = fields.index("channel")
        i_mpv = fields.index("mpv")
        i_error_mpv = fields.index("empv")

        for line in lines[2:]:
            v = line.split()
            mip_val = float(v[i_mpv])
            if float(v[i_error_mpv]) < 0:
                mip_val = 0
            idx_slab = self._N.slabs.index(int(v[i_slab]))
            mip_map[idx_slab][int(v[i_chip])][int(v[i_channel])] = mip_val
        if self._verbose:
            print("mip_map", mip_map)
        return mip_map


    def _read_masked(self, file_name):
        masked_map = np.zeros((
            self._N.n_slabs,
            self._N.n_chips,
            self._N.n_channels,
        ), dtype=int)
        print(aligned_path("Reading masked channels from ", file_name))
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
            idx_slab = self._N.slabs.index(int(v[0]))
            chip = int(v[1])
            for channel, mask_val in enumerate(v[2:]):
                masked_map[idx_slab][chip][channel] = mask_val
        if self._verbose:
            print("masked_map", masked_map)
        return masked_map


    def _handle_uncommissioned_positions(self, pedestal_map, mip_map, masked_map):
        """This changes the passed arrays in-place."""
        is_commissioned_map = np.ones_like(pedestal_map, dtype=int)
        is_commissioned_map[masked_map == 1] = 0

        # Handle pedestals
        sca_has_bad_pedestal = pedestal_map < self._N.pedestal_min_value
        if np.any(sca_has_bad_pedestal):
            sca_is_used_for_average = pedestal_map > self._N.pedestal_min_average
            channel_can_provide_average = np.expand_dims(
                sca_is_used_for_average.sum(axis=-1) >= self._N.pedestal_min_scas,
                sca_is_used_for_average.ndim - 1,
            )
            pedestal_has_no_fix = np.logical_and(
                sca_has_bad_pedestal,
                np.logical_not(channel_can_provide_average),
            )
            average_pedestal = np.empty_like(pedestal_map)
            average_pedestal[:] = np.expand_dims(
                np.divide(
                    np.multiply(pedestal_map, sca_is_used_for_average).sum(axis=-1),
                    np.maximum(1, sca_is_used_for_average.sum(axis=-1)),
                ),
                sca_is_used_for_average.ndim - 1,
            )
            pedestal_map[sca_has_bad_pedestal] = average_pedestal[sca_has_bad_pedestal]

            is_commissioned_map[np.logical_and(
                sca_has_bad_pedestal,
                channel_can_provide_average,
            )] = 0  # We might want to consider this case as "is_commissioned".
            is_commissioned_map[pedestal_has_no_fix] = 0
            pedestal_map[pedestal_has_no_fix] = -1
        if self._verbose:
            print("corrected pedestal_map", pedestal_map)

        # Handle MIPs
        has_bad_mip = mip_map < self._N.mip_cutoff
        if np.any(has_bad_mip):
            channel_is_used_for_average = mip_map >= self._N.mip_cutoff
            per_chip_average = channel_is_used_for_average.mean(axis=-1)
            per_chip_average[per_chip_average == 0] = self._N.mip_malfunctioning_chip
            mip_average_on_chip = np.empty_like(mip_map)
            mip_average_on_chip[:] = np.expand_dims(
                per_chip_average,
                per_chip_average.ndim,
            )
            mip_average_on_chip[np.isnan(mip_average_on_chip)] = -1
            is_commissioned_map[has_bad_mip] = 0
            mip_map[has_bad_mip] = mip_average_on_chip[has_bad_mip]
        if self._verbose:
            print("corrected mip_map", mip_map)

        if self._verbose:
            print("is_commissioned_map", is_commissioned_map)
            print("Rate of commissioned scas: {:.2f}%".format(is_commissioned_map.mean() * 100))
        return is_commissioned_map


def speed_warning_if_python2():
    if sys.version_info.major == 2:
        print(
            "Warning: Slow. The eventbuilding is run with python2. "
            "While the code aims to stay python2/3 compatible, "
            "it was found to run significantly faster with python3 "
            "(x5, see https://github.com/SiWECAL-TestBeam/SiWECAL-TB-analysis/pull/20#issuecomment-982763034)"
        )


if __name__ == "__main__":
    numbers = EcalNumbers(slabs=np.arange(15).tolist())
    EcalConfig(verbose=True, numbers=numbers)
