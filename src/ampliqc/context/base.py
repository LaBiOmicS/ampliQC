"""
Base context class and evaluation status definitions.
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Dict, Any, List


class Status(Enum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"


class QualityCheckResult:
    def __init__(
        self,
        name: str,
        status: Status,
        message: str,
        metrics: Dict[str, Any],
        context_reasoning: str,
    ):
        self.name = name
        self.status = status
        self.message = message
        self.metrics = metrics
        self.context_reasoning = context_reasoning

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status.value,
            "message": self.message,
            "metrics": self.metrics,
            "context_reasoning": self.context_reasoning,
        }


class ContextAnalyzer(ABC):
    """
    Abstract base class for context-aware sequencing data evaluations.
    """

    def __init__(self, context_name: str):
        self.context_name = context_name

    @abstractmethod
    def evaluate(self, stats: Dict[str, Any], sample_sequences: List[str]) -> Dict[str, Any]:
        """
        Runs context-specific evaluation on accumulated metrics and sample sequences.
        """
        pass
