"""Provide helper functions to execute command-line programs
concurrently."""

import asyncio
import dataclasses
import os
import textwrap
import threading
from typing import Any, Coroutine, Dict, List, Optional, TypeVar

from .helpers import PathLike, transient_dir

T = TypeVar("T")


@dataclasses.dataclass
class Command:
    """Represent a command-line program to be executed."""

    jobid: str
    """Job identifier. Important for tracking and error reporting if enabled."""
    args: List[str]
    """List of command-line arguments to execute the subprocess."""
    workdir: Optional[PathLike] = None
    """Working directory. If present, the subprocess will be run within
    the given directory. The directory will be created if not exists."""
    data: Dict[str, Any] = dataclasses.field(default_factory=dict)
    """Dictionary containing dynamic data that may be useful for
    handling the command."""


def async_run(coro: Coroutine[Any, Any, T]) -> T:
    """Execute the given coroutine in a run loop.

    Similar to :py:func:`asyncio.run` but reuses an existing run loop or
    creates a new one if needed, and returns the result of invoking the
    coroutine.

    .. note::
        This is a hack to run on a Jupyter Notebook, which includes its
        own run loop.
    """

    class RunThread(threading.Thread):
        def run(self):
            self.retvalue = asyncio.run(coro)

    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        thread = RunThread()
        thread.start()
        thread.join()
        return thread.retvalue
    else:
        return asyncio.run(coro)


async def concurrent_subprocess(
    commands: List[Command], tasks: int = 1, quiet: bool = False
) -> bool:
    """Run the given commands concurrently.

    Uses the `asyncio` module to run concurrent workers (controlled by
    *tasks*) as coroutines, which run the next available command from a
    FIFO queue.

    Note that the standard output of the subprocess is captured, and
    printed only when an error occurred and if *quiet* is False.

    :param List[Command] commands: List of commands to run.
    :param int tasks: Number of parallel tasks to run. Defaults to 1.
    :param bool quiet: Do not print progress to standard output.
        Defaults to False.
    :return bool: True if all subprocesses succeeded, else False.
    """

    async def worker():
        while True:
            cmd = await queue.get()
            i = n_cmd - queue.qsize()
            istr = "[{{:>{}}}/{{}}]".format(len(str(n_cmd))).format(i, n_cmd)
            if not quiet:
                print(f"{istr} Starting {cmd.jobid}...")
            with transient_dir(cmd.workdir or os.getcwd()):
                proc = await asyncio.create_subprocess_exec(
                    cmd.args[0],
                    *cmd.args[1:],
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
                stdout, _ = await proc.communicate()
                queue.task_done()
                if proc.returncode == 0:
                    continue

                await failed_queue.put(cmd)
                if not quiet:
                    output = textwrap.indent(
                        stdout.decode().strip(), " " * (len(istr) + 1)
                    )
                    print(f"{istr} Failed {cmd.jobid} with output:\n{output}")

    queue: asyncio.Queue[Command] = asyncio.Queue()
    failed_queue: asyncio.Queue[Command] = asyncio.Queue()
    n_cmd = len(commands)
    for cmd in commands:
        queue.put_nowait(cmd)

    if not quiet:
        print(f"Running {n_cmd} command(s) in {tasks} parallel task(s)...")
    workers: List[asyncio.Task[Any]] = [
        asyncio.create_task(worker()) for _ in range(tasks)
    ]
    await queue.join()

    for task in workers:
        task.cancel()
    await asyncio.gather(*workers, return_exceptions=True)

    success = failed_queue.qsize() == 0
    if not quiet and not success:
        print("WARNING: Some commands have failed. Please check the output.")
    return success
