"""Base execution interface for LLM operations."""

from abc import ABC, abstractmethod
from typing import List
from ..chunking.chunker import Chunk
from ..config import Config


class BaseExecutor(ABC):
    """Abstract base class for LLM execution modes."""
    
    def __init__(self, config: Config):
        self.config = config
    
    @abstractmethod
    def initialize(self) -> None:
        """Initialize the execution environment."""
        pass
    
    @abstractmethod
    def execute_completion(self, prompt: str) -> str:
        """Execute a completion request and return the result."""
        pass
    
    @abstractmethod
    def get_execution_info(self) -> dict:
        """Get information about the current execution configuration."""
        pass
