"""Coalescer that concatenates summaries."""

from typing import List
from datetime import datetime
from pathlib import Path

from ..config import Config
from ..utils.logger import setup_logger

logger = setup_logger(__name__)


class Coalescer:
    """Coalescer that concatenates summaries."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def coalesce_summaries(self, summaries: List[str], source_file: str = "unknown") -> str:
        """Combine summaries into final output."""
        logger.info(f"Coalescing {len(summaries)} summaries")
        
        # Create header
        header = self._create_header(source_file, len(summaries))
        
        # Join all summaries
        content = "\n\n---\n\n".join(summaries)
        
        # Combine header and content
        final_output = f"{header}\n\n{content}"
        
        logger.info(f"Created final summary: {len(final_output)} characters")
        return final_output
    
    def _create_header(self, source_file: str, num_chunks: int) -> str:
        """Create a header for the final summary."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Get the correct model name based on the active provider
        provider = self.config.llm_provider.lower()
        if provider == "groq":
            model_name = self.config.groq_model
        elif provider == "together":
            model_name = self.config.together_model
        elif provider == "anthropic":
            model_name = self.config.anthropic_model
        elif provider == "huggingface":
            model_name = self.config.hf_model_name
        else:
            model_name = f"{provider} (unknown model)"
        
        return f"""# Mathematics Textbook Summary

**Source:** {source_file}  
**Generated:** {timestamp}  
**Provider:** {self.config.llm_provider.upper()}  
**Model:** {model_name}  
**Chunks processed:** {num_chunks}  
**Tool:** Math Summarizer v1.0.0

---"""
    
    def save_summary(self, summary: str, output_path: Path) -> None:
        """Save the final summary to a file."""
        logger.info(f"Saving summary to: {output_path}")
        
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(summary)
            logger.info(f"Successfully saved summary to {output_path}")
        except Exception as e:
            logger.error(f"Error saving summary: {str(e)}")
            raise
