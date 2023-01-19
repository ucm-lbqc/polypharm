import contextlib
import os
import sys
from importlib import resources
from pathlib import Path
from typing import Any, Generator, Union

import jinja2

PathLike = Union[str, Path]

ROOT_DIR = Path(__file__).parent.parent
TEMPLATE_ENV = jinja2.Environment(
    loader=jinja2.FileSystemLoader(ROOT_DIR / "polypharm" / "templates")
)


def get_schrodinger_path() -> str:
    """Return the ``SCHRODINGER_PATH`` environment variable.

    :raises OSError: ``SCHRODINGER_PATH`` is not set.
    :return str: The path in ``SCHRODINGER_PATH``.
    """
    path = os.getenv("SCHRODINGER_PATH")
    if not path:
        raise OSError("Environment variable SCHRODINGER_PATH is not set")
    if not os.path.exists(path):
        print("error: SCHRODINGER_PATH is invalid")
    return path


def get_script_path(filename: str) -> str:
    """Return the path to the given script.

    Scripts are located in the ``scripts`` directory within the package.

    :param str filename: Filename of the script (including extension).
    :raises FileNotFoundError: Script does not exist.
    :return str: Absolute path to the script.
    """
    resource = resources.files("polypharm").joinpath(os.path.join("scripts", filename))
    with resources.as_file(resource) as script_file:
        if not script_file.exists():
            raise FileNotFoundError(f"Script {filename} not found")
        return str(script_file)


def render_template(filename: str, **vars: Any) -> str:
    """Render the given template to a string.

    A template is a plain text file containing placeholders to be
    substituted with dynamic data on runtime. Uses the Jinja_ template
    engine for the rendering.

    :param str filename: Filename of the template (including extension).
    :param \\**vars: Dictionary containing the variables and the
        corresponding values to replace the placeholders in the
        template.
    :return str: The rendered template with the given data.

    .. _Jinja: https://jinja.palletsprojects.com/
    """
    template = TEMPLATE_ENV.get_template(filename)
    return template.render(vars)


def render_template_to_file(template_file: str, path: PathLike, **vars: Any) -> None:
    """Convenience function to render the given template to a file.

    :param str filename: Filename of the template (including extension).
    :param str | Path path: Path to the file to write the rendered
        template to.
    :param \\**vars: Dictionary containing the variables and the
        corresponding values to replace the placeholders in the
        template.
    """
    with open(path, "w") as io:
        io.write(render_template(template_file, **vars))


@contextlib.contextmanager
def transient_dir(workdir: PathLike) -> Generator[None, None, None]:
    """Context manager that executes code in the given directory, and
    then restore the original directory.

    If the working directory or any of its parents does not exists, they
    will be created.

    :param str | Path workdir: Working directory to change to.
    """
    Path(workdir).mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(workdir)
        yield
    finally:
        os.chdir(cwd)
