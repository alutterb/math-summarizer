"""Tests for the chunking module."""

import pytest
from pathlib import Path

from src.math_summarizer.config import Config
from src.math_summarizer.chunking.markdown_chunker import MarkdownChunker
from src.math_summarizer.utils.math_parser import MathParser


class TestMarkdownChunker:
    """Test the MarkdownChunker class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.config = Config(
            openai_api_key="test-key",
            max_chunk_size=1000,
            chunk_overlap=100
        )
        self.chunker = MarkdownChunker(self.config)
    
    def test_simple_chunking(self):
        """Test basic chunking functionality."""
        content = """
# Chapter 1: Introduction

This is a simple introduction to mathematics.

## 1.1 Basic Concepts

Here we define basic mathematical concepts.

### Definition 1.1.1
A **set** is a collection of objects.

## 1.2 Examples

Let's look at some examples.
"""
        chunks = self.chunker.chunk_content(content, "test.md")
        
        assert len(chunks) >= 1
        assert all(chunk.token_count <= self.config.max_chunk_size for chunk in chunks)
    
    def test_math_detection(self):
        """Test detection of mathematical content."""
        content = """
# Mathematical Content

Here is some inline math: $x + y = z$.

And here is display math:
$$\\int_0^1 f(x) dx = F(1) - F(0)$$

### Theorem 1.1
For all real numbers $a$ and $b$, we have $a + b = b + a$.

**Proof**: This follows from the commutative property of addition.
"""
        chunks = self.chunker.chunk_content(content, "test.md")
        
        # Check that mathematical content is detected
        math_chunks = [chunk for chunk in chunks if chunk.contains_math]
        assert len(math_chunks) > 0
        
        theorem_chunks = [chunk for chunk in chunks if chunk.contains_theorem]
        proof_chunks = [chunk for chunk in chunks if chunk.contains_proof]
        
        # Should detect theorem and proof content
        assert len(theorem_chunks) > 0 or len(proof_chunks) > 0


class TestMathParser:
    """Test the MathParser class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.parser = MathParser()
    
    def test_inline_math_detection(self):
        """Test detection of inline math."""
        text = "Here is some math: $x^2 + y^2 = z^2$ in the text."
        elements = self.parser.parse_mathematical_elements(text)
        
        inline_math = [e for e in elements if e.element_type == 'inline_math']
        assert len(inline_math) == 1
        assert inline_math[0].metadata['math_content'] == 'x^2 + y^2 = z^2'
    
    def test_display_math_detection(self):
        """Test detection of display math."""
        text = """
Here is display math:
$$\\sum_{i=1}^n i = \\frac{n(n+1)}{2}$$
End of math.
"""
        elements = self.parser.parse_mathematical_elements(text)
        
        display_math = [e for e in elements if e.element_type == 'display_math']
        assert len(display_math) == 1
        assert 'sum' in display_math[0].metadata['math_content']
    
    def test_theorem_detection(self):
        """Test detection of theorem structures."""
        text = """
### Theorem 1.2.3 (Fundamental Theorem)
Every continuous function on a closed interval is uniformly continuous.
"""
        structure = self.parser.identify_theorem_structure(text)
        
        assert structure is not None
        assert structure['type'] == 'theorem'
        assert structure['number'] == '1.2.3'
        assert structure['name'] == 'Fundamental Theorem'
    
    def test_proof_detection(self):
        """Test detection of proof structures."""
        text = """
**Proof**: We proceed by contradiction. Assume that...
Therefore, the statement holds. ∎
"""
        structure = self.parser.identify_proof_structure(text)
        
        assert structure is not None
        assert 'contradiction' in structure['content'].lower()
    
    def test_safe_split_points(self):
        """Test identification of safe split points."""
        text = "Before math $x + y$ after math. More text here."
        
        # Position inside math should not be safe
        assert not self.parser.is_safe_split_point(text, 15)  # Inside $x + y$
        
        # Position outside math should be safe
        assert self.parser.is_safe_split_point(text, 25)  # After math
