"""Anthropic Claude LLM client implementation."""

from .base import BaseLLMClient
import logging

logger = logging.getLogger(__name__)


class AnthropicClient(BaseLLMClient):
    """Anthropic Claude client implementation."""
    
    def _initialize_client(self) -> None:
        """Initialize the Anthropic client."""
        try:
            import anthropic
            
            
            self._client = anthropic.Anthropic(
                api_key=self.api_key,
                timeout=self.timeout
            )
            logger.info(f"Anthropic client initialized with model: {self.model}")
            
        except ImportError:
            raise ImportError("anthropic is required for Anthropic support. Install with: pip install anthropic")
        except Exception as e:
            logger.error(f"Failed to initialize Anthropic client: {str(e)}")
            raise
    
    def _generate_text(self, prompt: str) -> str:
        """Generate text using Anthropic API."""
        try:
            response = self._client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                temperature=self.temperature,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )
            return response.content[0].text
        except Exception as e:
            logger.error(f"Error generating text with Anthropic: {str(e)}")
            raise
