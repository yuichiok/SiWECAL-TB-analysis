import sys


class EventBuildingException(Exception):
    pass


def aligned_path(text, path):
    _text_length = 29
    assert _text_length >= len(text)
    return text + (_text_length - len(text)) * " " + str(path)


def speed_warning_if_python2():
    if sys.version_info.major == 2:
        print(
            "Warning: Slow. The eventbuilding is run with python2. "
            "While the code aims to stay python2/3 compatible, "
            "it was found to run significantly faster with python3 "
            "(x5, see https://github.com/SiWECAL-TestBeam/SiWECAL-TB-analysis/pull/20#issuecomment-982763034)"
        )


try:
    from tqdm.autonotebook import tqdm

    def get_tree_spills(tree, max_entries):
        for i, spill in enumerate(tqdm(tree, desc="# Build events", total=max_entries, unit=" spills")):
            if i > max_entries:
                break
            yield i, spill
except ImportError:
    from datetime import datetime

    def get_tree_spills(tree, max_entries):
        print("# Going to analyze %i entries..." %max_entries)
        print("# For better progress information: `pip install tqdm`.")
        print("# Start time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        progress_bar = ""
        for i, spill in enumerate(tree):
            if i > max_entries:
                break
            if i%10 == 0 or i == max_entries:
                progress_bar = "#" * int(30 * i / max_entries)
                print("# Build events [{}] Spill {}/{}".format(
                        progress_bar.ljust(30),
                        str(i).rjust(len(str(max_entries))),
                        max_entries,
                    ),
                    end="\r",
                )
                sys.stdout.flush()
            yield i, spill
        print("\n# Final time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def get_tree_spills_no_progress_info(tree, max_entries):
    for i, spill in enumerate(tree):
        if i > max_entries:
            break
        yield i, spill
