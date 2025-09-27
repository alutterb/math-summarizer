"""API-based execution for remote LLM services."""

import logging
from .base import BaseExecutor
from ..llm_clients import LLMClientFactory, BaseLLMClient
from ..config import Config

logger = logging.getLogger(__name__)


class APIExecutor(BaseExecutor):
    """Handles execution using remote API services."""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.main_client: BaseLLMClient = None
    
    def initialize(self) -> None:
        """Initialize API clients."""
        logger.info(f"Initializing API mode with provider: {self.config.llm_provider}")
        
        try:
            # Create main client
            self.main_client = LLMClientFactory.create_main_client(self.config)
            logger.info(f"Main client initialized: {self.main_client.get_config_summary()}")
            
        except Exception as e:
            logger.error(f"Failed to initialize API clients: {str(e)}")
            raise
    
    def execute_completion(self, prompt: str) -> str:
        """Execute completion using the API client."""
        if self.main_client is None:
            raise RuntimeError("API executor not initialized. Call initialize() first.")
        
        return self.main_client.complete(prompt)
    
    def get_execution_info(self) -> dict:
        """Get API execution information."""
        if self.main_client is None:
            return {"mode": "api", "status": "not_initialized"}
        
        return {
            "mode": "api",
            "provider": self.config.llm_provider,
            "model": self.main_client.model,
            "status": "initialized"
        }
