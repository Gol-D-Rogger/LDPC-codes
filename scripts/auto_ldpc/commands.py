"""
Command execution helpers with optional WSL prefix and dry-run support.
"""

from __future__ import annotations

import shlex
import subprocess
from typing import Iterable, List, Optional, Sequence


class CommandRunner:
    def __init__(
        self,
        command_prefix: Optional[Sequence[str]] = None,
        dry_run: bool = False,
    ) -> None:
        self.command_prefix = list(command_prefix or [])
        self.dry_run = dry_run

    def build_command(self, args: Sequence[str]) -> List[str]:
        return self.command_prefix + list(args)

    def run(
        self,
        args: Sequence[str],
        capture_output: bool = False,
        check: bool = True,
        text: bool = True,
        cwd: Optional[str] = None,
        env: Optional[dict] = None,
    ) -> subprocess.CompletedProcess:
        cmd = self.build_command(args)
        if self.dry_run:
            printable = " ".join(shlex.quote(part) for part in cmd)
            print(f"[DRY-RUN] {printable}")
            return subprocess.CompletedProcess(cmd, 0, "", "")
        return subprocess.run(
            cmd,
            capture_output=capture_output,
            text=text,
            check=check,
            cwd=cwd,
            env=env,
        )
