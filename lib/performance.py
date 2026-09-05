"""Wall-time accounting for the tasks a production run performs.

OpenMM reporters are duck-typed: anything exposing ``describeNextReport`` and
``report`` is one. Timing therefore needs no change to OpenMM itself. Each real
reporter is wrapped in a proxy that times its ``report`` call, the workflow's
own checkpoint and monitor calls are timed directly, and whatever is left over
is integration.
"""

from __future__ import annotations

from contextlib import contextmanager
import csv
from pathlib import Path
import time


PERFORMANCE_FIELDS = (
    "production_time_ns",
    "step",
    "wall_s",
    "interval_s",
    "ns_per_day",
    "md_s",
    "trajectory_s",
    "state_s",
    "checkpoint_s",
    "monitor_s",
    "reported_pct",
)
TASKS = ("trajectory", "state", "checkpoint", "monitor")


class PerformanceTracker:
    """Accumulate wall time per task since production began.

    Timings cover the current submission only: a restart is a new process, so
    its accounting starts from zero again.
    """

    def __init__(self) -> None:
        self.totals: dict[str, float] = {task: 0.0 for task in TASKS}
        self._started: float | None = None

    def start(self) -> None:
        self._started = time.perf_counter()

    @property
    def started(self) -> bool:
        return self._started is not None

    def elapsed(self) -> float:
        return 0.0 if self._started is None else time.perf_counter() - self._started

    @contextmanager
    def measure(self, task: str):
        """Attribute the enclosed work to one task."""
        start = time.perf_counter()
        try:
            yield
        finally:
            self.totals[task] = self.totals.get(task, 0.0) + (
                time.perf_counter() - start
            )

    def summary(self) -> str:
        """Render the final breakdown for the end of a run."""
        wall = self.elapsed()
        reported = sum(self.totals.values())
        lines = [
            f"Production wall time: {wall:.1f} s "
            f"({reported / wall * 100:.2f}% in reporting and checkpointing)"
            if wall > 0
            else "Production wall time: 0.0 s"
        ]
        for task in TASKS:
            seconds = self.totals.get(task, 0.0)
            share = seconds / wall * 100 if wall > 0 else 0.0
            lines.append(f"  {task:<12} {seconds:8.2f} s  {share:6.2f}%")
        integration = wall - reported
        share = integration / wall * 100 if wall > 0 else 0.0
        lines.append(f"  {'integration':<12} {integration:8.2f} s  {share:6.2f}%")
        return "\n".join(lines)


class TimedReporter:
    """Delegate to a real reporter while timing what its report costs.

    OpenMM builds one ``State`` for every reporter due at the same step before
    any of them is called, so that shared cost stays in integration rather than
    being attributed here.
    """

    def __init__(self, reporter: object, tracker: PerformanceTracker, task: str) -> None:
        self._reporter = reporter
        self._tracker = tracker
        self._task = task

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM's interface
        return self._reporter.describeNextReport(simulation)

    def report(self, simulation, state):
        with self._tracker.measure(self._task):
            self._reporter.report(simulation, state)

    def __getattr__(self, name):
        return getattr(self._reporter, name)


class PerformanceReporter:
    """Write one timing row every ``interval_steps`` of production."""

    def __init__(
        self,
        path: Path,
        interval_steps: int,
        tracker: PerformanceTracker,
        timestep_ns: float,
        append: bool = False,
        initial_step: int = 0,
    ) -> None:
        self._path = Path(path)
        self._interval = interval_steps
        self._tracker = tracker
        self._timestep_ns = timestep_ns
        # Rates are measured against where this submission actually started, so
        # the first row after a restart is not computed from time zero.
        self._previous = (initial_step * timestep_ns, 0.0)
        if not append or not self._path.is_file() or self._path.stat().st_size == 0:
            with self._path.open("w", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=PERFORMANCE_FIELDS).writeheader()

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM's interface
        steps = self._interval - simulation.currentStep % self._interval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        wall = self._tracker.elapsed()
        production_ns = simulation.currentStep * self._timestep_ns
        previous_ns, previous_wall = self._previous
        interval_s = wall - previous_wall
        ns_per_day = (
            (production_ns - previous_ns) / interval_s * 86400 if interval_s > 0 else 0.0
        )
        self._previous = (production_ns, wall)
        reported = sum(self._tracker.totals.values())
        row = {
            "production_time_ns": f"{production_ns:.12g}",
            "step": simulation.currentStep,
            "wall_s": f"{wall:.3f}",
            "interval_s": f"{interval_s:.3f}",
            "ns_per_day": f"{ns_per_day:.3f}",
            "md_s": f"{wall - reported:.3f}",
            "trajectory_s": f"{self._tracker.totals.get('trajectory', 0.0):.3f}",
            "state_s": f"{self._tracker.totals.get('state', 0.0):.3f}",
            "checkpoint_s": f"{self._tracker.totals.get('checkpoint', 0.0):.3f}",
            "monitor_s": f"{self._tracker.totals.get('monitor', 0.0):.3f}",
            "reported_pct": f"{reported / wall * 100:.3f}" if wall > 0 else "0",
        }
        with self._path.open("a", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=PERFORMANCE_FIELDS).writerow(row)
