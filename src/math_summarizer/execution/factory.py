"""Factory for creating execution mode handlers."""

import logging
from .base import BaseExecutor
from .api_executor import APIExecutor
from .local_executor import LocalExecutor
from ..config import Config

logger = logging.getLogger(__name__)


class ExecutionFactory:
    """Factory for creating execution mode handlers."""
    
    @staticmethod
    def create_executor(config: Config) -> BaseExecutor:
        """Create the appropriate executor based on configuration."""
        mode = config.execution_mode.lower()
        
        if mode == "api":
            executor = APIExecutor(config)
        elif mode == "local":
            executor = LocalExecutor(config)
        else:
            raise ValueError(f"Unsupported execution mode: {mode}. Supported modes: api, local")
        
        # Initialize the executor
        executor.initialize()
        
        logger.info(f"Created and initialized {mode} executor: {executor.get_execution_info()}")
        return executor
