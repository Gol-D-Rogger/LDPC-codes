from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class IntervalPlan:
    start: float
    stop: float
    axis_type: str
    direction: int  # +1 means increasing axis, -1 means decreasing axis


@dataclass(frozen=True)
class StepPolicy:
    step_default: float
    step_mid: float
    step_low: float
    fer_mid: float
    fer_low: float

    def step_for(self, fer: Optional[float]) -> float:
        if fer is None or not math.isfinite(fer):
            return self.step_default
        if fer <= self.fer_low:
            return self.step_low
        if fer <= self.fer_mid:
            return self.step_mid
        return self.step_default
