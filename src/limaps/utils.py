import logging

logger = logging.getLogger(__name__)
import json
from pathlib import Path

try:
    import tkinter
    import tkinter.filedialog as dialog

except ImportError as e:
    logger.warn(e)


def askdirectory() -> Path:
    root = tkinter.Tk()
    root.focus_force()
    root.wm_withdraw()
    targetfile = (
        dialog.askopenfilename(
            parent=root,
            initialdir=Path("."),
            filetypes=[
                ("Format text(.csv, .xls, .xlsx)", "*.csv *.xls *xlsx"),
                ("Plain text", "*.txt"),
            ],
        )
        or ""
    )
    if not targetfile:
        return None
    return Path(targetfile)


def save_figure(figpath, fig, transparent=True, dpi=150, **kwargs):
    figpath = Path(figpath)
    if not callable(getattr(fig, "savefig", None)):
        err = "Input is not a matplotlib.figure.Figure instance"
        raise TypeError(err)
    fig.savefig(figpath, transparent=transparent, dpi=dpi, **kwargs)
    logger.info(f"Figure saved at {figpath}.")


def load_json(path):
    with open(path, mode="r") as file:
        data = json.load(file)
    return data
