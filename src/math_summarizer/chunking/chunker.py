"""Chunker that splits content by lines or tokens."""

from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass

from ..config import Config
from ..utils.logger import setup_logger

logger = setup_logger(__name__)


@dataclass
class Chunk:
    """Represents a chunk of content."""
    content: str
    chunk_id: str
    chunk_index: int
    token_count: int
    source: str


class Chunker:
    """Chunker that splits content by lines or tokens."""
    
    def __init__(self, config: Config):
        self.config = config
        
        # Initialize tokenizer for the Hugging Face model
        try:
            if config.execution_mode.lower() == "api":
                # For API mode, try to use tokenizers library if available
                try:
                    from tokenizers import Tokenizer
                    from huggingface_hub import hf_hub_download
                    
                    # Try to download tokenizer
                    tokenizer_path = hf_hub_download(
                        repo_id=config.hf_model_name,
                        filename="tokenizer.json",
                        token=config.hf_token
                    )
                    self.tokenizer = Tokenizer.from_file(tokenizer_path)
                    logger.info(f"Initialized tokenizer for {config.hf_model_name} (API mode)")
                except Exception as e:
                    logger.info(f"Could not load tokenizer in API mode: {str(e)}, using fallback")
                    self.tokenizer = None
            else:
                # Local mode - use transformers
                from transformers import AutoTokenizer
                self.tokenizer = AutoTokenizer.from_pretrained(
                    config.hf_model_name,
                    token=config.hf_token,
                    trust_remote_code=True
                )
                logger.info(f"Initialized tokenizer for {config.hf_model_name} (local mode)")
        except Exception as e:
            logger.warning(f"Could not load tokenizer: {str(e)}")
            logger.info("Using fallback token counting (approximate)")
            self.tokenizer = None
    
    def chunk_file(self, file_path: Path) -> List[Chunk]:
        """Chunk a file into simple chunks."""
        logger.info(f"Chunking file: {file_path}")
        
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        return self.chunk_content(content, str(file_path))
    
    def chunk_content(self, content: str, source: str = "unknown") -> List[Chunk]:
        """Split content into chunks based on token count."""
        chunks = []
        
        # Split by paragraphs first to avoid breaking mid-sentence
        paragraphs = content.split('\n\n')
        
        current_chunk = ""
        chunk_index = 0
        
        for paragraph in paragraphs:
            # Check if adding this paragraph would exceed token limit
            potential_chunk = current_chunk + "\n\n" + paragraph if current_chunk else paragraph
            token_count = self._count_tokens(potential_chunk)
            
            if token_count > self.config.max_chunk_size and current_chunk:
                # Save current chunk and start new one
                chunk = self._create_chunk(current_chunk, source, chunk_index)
                chunks.append(chunk)
                chunk_index += 1
                current_chunk = paragraph
            else:
                current_chunk = potential_chunk
        
        # Add final chunk
        if current_chunk.strip():
            chunk = self._create_chunk(current_chunk, source, chunk_index)
            chunks.append(chunk)
        
        logger.info(f"Created {len(chunks)} chunks from {source}")
        return chunks
    
    def chunk_by_lines(self, content: str, source: str = "unknown", lines_per_chunk: int = 100) -> List[Chunk]:
        """Split content by number of lines."""
        lines = content.split('\n')
        chunks = []
        
        for i in range(0, len(lines), lines_per_chunk):
            chunk_lines = lines[i:i + lines_per_chunk]
            chunk_content = '\n'.join(chunk_lines)
            
            if chunk_content.strip():  # Skip empty chunks
                chunk = self._create_chunk(chunk_content, source, len(chunks))
                chunks.append(chunk)
        
        logger.info(f"Created {len(chunks)} chunks from {source} (by lines)")
        return chunks
    
    def _create_chunk(self, content: str, source: str, chunk_index: int) -> Chunk:
        """Create a Chunk."""
        token_count = self._count_tokens(content)
        chunk_id = f"{source}_chunk_{chunk_index}"
        
        return Chunk(
            content=content.strip(),
            chunk_id=chunk_id,
            chunk_index=chunk_index,
            token_count=token_count,
            source=source
        )
    
    def _count_tokens(self, text: str) -> int:
        """Count tokens in text using the model's tokenizer."""
        if self.tokenizer is not None:
            try:
                if hasattr(self.tokenizer, 'encode'):
                    # Transformers tokenizer
                    if hasattr(self.tokenizer, 'add_special_tokens'):
                        return len(self.tokenizer.encode(text, add_special_tokens=False))
                    else:
                        return len(self.tokenizer.encode(text))
                elif hasattr(self.tokenizer, 'encode'):
                    # Tokenizers library tokenizer
                    return len(self.tokenizer.encode(text).ids)
                else:
                    logger.warning("Unknown tokenizer type, using fallback")
            except Exception as e:
                logger.warning(f"Error counting tokens: {str(e)}, using fallback")
        
        # Fallback: approximate token count (roughly 4 characters per token)
        return len(text) // 4
