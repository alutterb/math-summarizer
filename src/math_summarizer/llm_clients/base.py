"""Base abstract class for LLM clients."""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import time
import logging

logger = logging.getLogger(__name__)


class BaseLLMClient(ABC):
    """Abstract base class for all LLM clients."""
    
    def __init__(self, model: str, api_key: str, temperature: float = 0.7, 
                 max_tokens: int = 512, timeout: int = 120, max_retries: int = 3):
        self.model = model
        self.api_key = api_key
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.timeout = timeout
        self.max_retries = max_retries
        self._client = None
        
    @abstractmethod
    def _initialize_client(self) -> None:
        """Initialize the provider-specific client."""
        pass
    
    @abstractmethod
    def _generate_text(self, prompt: str) -> str:
        """Generate text using the provider's API."""
        pass
    
    def complete(self, prompt: str) -> str:
        """Complete text generation with retry logic."""
        if self._client is None:
            self._initialize_client()
            
        for attempt in range(self.max_retries):
            try:
                logger.debug(f"Text generation attempt {attempt + 1}/{self.max_retries}")
                result = self._generate_text(prompt)
                return result.strip()
                
            except Exception as e:
                logger.warning(f"Generation attempt {attempt + 1} failed: {str(e)}")
                if attempt < self.max_retries - 1:
                    wait_time = 2 ** attempt
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"All {self.max_retries} attempts failed")
                    raise
    
    @property
    def provider_name(self) -> str:
        """Get the provider name."""
        return self.__class__.__name__.replace('Client', '').lower()
    
    def get_config_summary(self) -> Dict[str, Any]:
        """Get configuration summary for logging."""
        return {
            'provider': self.provider_name,
            'model': self.model,
            'temperature': self.temperature,
            'max_tokens': self.max_tokens,
            'timeout': self.timeout
        }
