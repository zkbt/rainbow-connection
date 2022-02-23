import matplotlib.animation as ani


def get_writer(filename, fps=30, **kw):
    """
    Try to get an appropriate animation writer,
    given the filename provided.

    Parameters
    ----------

    filename : str
            The output filename string for the animation.

    fps : float
            Frames/second.

    kw : dict
            All other keywords will be passed to the initialization
            of the animation writer.
    """
    if ".mp4" in filename:
        try:
            writer = ani.writers["ffmpeg"](fps=fps, **kw)
        except (RuntimeError, KeyError):
            message = """
            This computer seems unable to run ffmpeg
            to create `.mp4` movies. Try running
                `conda install ffmpeg`
            and trying again.
            """
            raise RuntimeError()
    else:
        try:
            writer = ani.writers["pillow"](fps=fps, **kw)
        except (RuntimeError, KeyError):
            writer = ani.writers["imagemagick"](fps=fps, **kw)
            raise RuntimeError("This computer seem unable to animate?")
    return writer
