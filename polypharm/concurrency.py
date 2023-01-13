import asyncio
import dataclasses
import os
import subprocess
import threading
from typing import Any, Coroutine, Dict, List, Optional

from .helpers import PathLike, transient_dir


@dataclasses.dataclass
class Command:
    jobid: str
    args: List[str]
    workdir: Optional[PathLike] = None
    data: Dict[str, Any] = dataclasses.field(default_factory=dict)


# hack to run on a Jupyter Notebook (use existing run loop)
def async_run(coro: Coroutine[Any, Any, None]) -> None:
    class RunThread(threading.Thread):
        def run(self):
            asyncio.run(coro)

    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        thread = RunThread()
        thread.start()
        thread.join()
    else:
        return asyncio.run(coro)


async def concurrent_subprocess(
    commands: List[Command], tasks: int = 1, quiet: bool = False
) -> None:
    async def worker():
        while True:
            cmd = await queue.get()
            if not quiet:
                i = nCommands - queue.qsize()
                print(f"[{i}/{nCommands}] Starting {cmd.jobid}...")
            with transient_dir(cmd.workdir or os.getcwd()):
                proc = await asyncio.create_subprocess_exec(
                    cmd.args[0],
                    *cmd.args[1:],
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
                stdout, stderr = await proc.communicate()
                if proc.returncode != 0:
                    raise subprocess.CalledProcessError(
                        proc.returncode or 0, cmd.args, stdout, stderr
                    )
                queue.task_done()

    queue: asyncio.Queue[Command] = asyncio.Queue()
    nCommands = len(commands)
    for cmd in commands:
        queue.put_nowait(cmd)

    if not quiet:
        print(f"Running {nCommands} command(s) in {tasks} parallel task(s)...")
    workers: List[asyncio.Task[Any]] = [
        asyncio.create_task(worker()) for _ in range(tasks)
    ]
    await queue.join()

    for task in workers:
        task.cancel()
    await asyncio.gather(*workers, return_exceptions=True)
