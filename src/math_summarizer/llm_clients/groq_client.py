"""Groq LLM client implementation."""

from .base import BaseLLMClient
import logging

logger = logging.getLogger(__name__)


class GroqClient(BaseLLMClient):
    """Groq API client implementation."""
    
    def _initialize_client(self) -> None:
        """Initialize the Groq client."""
        try:
            from llama_index.llms.groq import Groq
            
            self._client = Groq(
                model=self.model,
                api_key=self.api_key,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
            )
            logger.info(f"Groq client initialized with model: {self.model}")
            
        except ImportError:
            raise ImportError("llama-index-llms-groq is required for Groq support. Install with: pip install llama-index-llms-groq")
        except Exception as e:
            logger.error(f"Failed to initialize Groq client: {str(e)}")
            raise
    
    def _generate_text(self, prompt: str) -> str:
        """Generate text using Groq API."""
        response = self._client.complete(prompt)
        return response.text
