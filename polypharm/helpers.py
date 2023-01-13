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
    loader=jinja2.FileSystemLoader(ROOT_DIR / "templates")
)


def get_schrodinger_path() -> str:
    path = os.getenv("SCHRODINGER_PATH")
    if not path:
        print(
            "error: Environment variable SCHRODINGER_PATH is not set", file=sys.stderr
        )
        sys.exit(1)
    if not os.path.exists(path):
        print("error: SCHRODINGER_PATH is invalid")
    return path


def get_script_path(filename: str) -> str:
    resource = resources.files("polypharm")
    with resources.as_file(resource) as package_dir:
        script_file = package_dir.parent / "scripts" / filename
        if not script_file.exists():
            raise FileNotFoundError(f"Script {filename} not found")
        return str(script_file)


def render_template(filename: str, **vars: Any) -> str:
    template = TEMPLATE_ENV.get_template(filename)
    return template.render(vars)


def render_template_to_file(template_file: str, path: PathLike, **vars: Any) -> None:
    with open(path, "w") as io:
        io.write(render_template(template_file, **vars))


@contextlib.contextmanager
def transient_dir(path: PathLike) -> Generator[None, None, None]:
    Path(path).mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)
