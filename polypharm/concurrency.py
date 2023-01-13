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
    jobid: str
    args: List[str]
    workdir: Optional[PathLike] = None
    data: Dict[str, Any] = dataclasses.field(default_factory=dict)


# hack to run on a Jupyter Notebook (use existing run loop)
def async_run(coro: Coroutine[Any, Any, T]) -> T:
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
