"""Anthropic Claude LLM client implementation."""

from .base import BaseLLMClient
import logging

logger = logging.getLogger(__name__)


class AnthropicClient(BaseLLMClient):
    """Anthropic Claude client implementation."""
    
    def _initialize_client(self) -> None:
        """Initialize the Anthropic client."""
        try:
            from llama_index.llms.anthropic import Anthropic
            
            self._client = Anthropic(
                model=self.model,
                api_key=self.api_key,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
            )
            logger.info(f"Anthropic client initialized with model: {self.model}")
            
        except ImportError:
            raise ImportError("llama-index-llms-anthropic is required for Anthropic support. Install with: pip install llama-index-llms-anthropic")
        except Exception as e:
            logger.error(f"Failed to initialize Anthropic client: {str(e)}")
            raise
    
    def _generate_text(self, prompt: str) -> str:
        """Generate text using Anthropic API."""
        response = self._client.complete(prompt)
        return response.text
