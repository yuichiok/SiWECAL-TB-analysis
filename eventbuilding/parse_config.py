#!/usr/bin/env python
from __future__ import print_function

import argparse
import configparser
import os

_default_cfg_dir = os.path.dirname(os.path.abspath(__file__))
_default_cfg = os.path.join(_default_cfg_dir, "default_eventbuilding.cfg")
assert os.path.exists(_default_cfg), _default_cfg


def parse_config(config_file=None):
    if config_file is None:
        config_file = _default_cfg
    parser = configparser.ConfigParser()
    parser.read(config_file)
    return parser


def get_comments_from_config(config_file=None):
    if config_file is None:
        config_file = _default_cfg
    with open(config_file) as f:
        lines = f.readlines()
    comments = {None: {}}
    current_section = comments[None]
    current_comment = " "
    for line in map(str.strip, lines):
        if len(line.strip()) == 0:
            current_comment = " "
        elif line.startswith("["):
            assert line.endswith("]")
            current_section_name = line[1:-1]
            if current_section_name not in comments:
                comments[current_section_name] = {}
            current_section = comments[current_section_name]
        elif line.startswith("#") or line.startswith(";"):
            current_comment += line[2:]
        elif line.find("=") >= 0 or line.find(":") >= 0:
            find1 = line.find("=")
            find2 = line.find(":")
            if find1 >= 0 and find2 >= 0:
                find = min(find1, find2)
            else:
                find = max(find1, find2)
            current_section[line[:find].strip()] = current_comment
            current_comment = " "
        else:
            raise NotImplementedError(line, len(line))
    return comments


def _read_config_file(config_parser):
    base_arg_parser = argparse.ArgumentParser(add_help=False)
    _help = "Values specified here overwrite the default configuration. "
    _help += "CLI values take precedence over any configuration file settings"
    base_arg_parser.add_argument("-c", "--config_file", default=None, help=_help)
    args, remaining_argv = base_arg_parser.parse_known_args()
    if args.config_file:
        assert os.path.exists(args.config_file), args.config_file
        config_parser.read(args.config_file)
    return base_arg_parser, remaining_argv


def create_cli_from_default_config(config_file=None, section="eventbuilding"):
    config_parser = parse_config(config_file)
    base_arg_parser, remaining_argv = _read_config_file(config_parser)
    arg_parser = argparse.ArgumentParser(
        description="Build an event-level rootfile (smaller) from the raw rootfile.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[base_arg_parser],
    )
    comments = get_comments_from_config(config_file)
    for key, default in config_parser[section].items():
        _help = comments.get(section, {}).get(key, " ")
        arg_parser.add_argument("--" + key, default=default, help=_help)
    ap_args = arg_parser.parse_args(remaining_argv)
    for field, value in vars(ap_args).items():
        if field == "config_file":
            continue
        config_parser[section][field] = value
    return config_parser


if __name__ == "__main__":
    config_parser = create_cli_from_default_config()
    for section in config_parser.sections():
        print("SECTION:", section)
        for k, v in config_parser[section].items():
            print(" -", k, v)
