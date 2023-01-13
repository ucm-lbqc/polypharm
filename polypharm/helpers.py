import contextlib
import os
import sys
from pathlib import Path
from typing import Any, Generator, Union

import jinja2

PathLike = Union[str, Path]

ROOT_DIR = Path(__file__).parent.parent
SCRIPT_DIR = ROOT_DIR / "scripts"
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
    path = Path(SCRIPT_DIR, filename)
    if not path.exists():
        raise FileNotFoundError(f"Script {filename} not found")
    return str(path)


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
