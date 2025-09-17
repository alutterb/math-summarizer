"""Simplified main CLI for the Math Summarizer."""

import click
from pathlib import Path
from typing import Optional
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .config import config
from .utils.logger import setup_logger
# Model utils not needed for HF models
from .chunking.chunker import Chunker
from .summarization.summarizer import Summarizer
from .coalescing.coalescer import Coalescer

console = Console()
logger = setup_logger(__name__)


@click.command()
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option('--output', '-o', type=click.Path(path_type=Path), 
              help='Output file path (default: input_file_summary.md)')
@click.option('--chunk-size', type=int, default=None,
              help=f'Maximum chunk size in tokens (default: {config.max_chunk_size})')
@click.option('--chunk-method', type=click.Choice(['tokens', 'lines']), default='tokens',
              help='Chunking method: by tokens or by lines (default: tokens)')
@click.option('--lines-per-chunk', type=int, default=100,
              help='Lines per chunk when using line-based chunking (default: 100)')
@click.option('--provider', type=click.Choice(['groq', 'together', 'anthropic', 'huggingface']), default=None,
              help=f'LLM provider to use (default: {config.llm_provider})')
@click.option('--model', type=str, default=None,
              help='Model name to use (provider-specific)')
@click.option('--temperature', type=float, default=None,
              help=f'Temperature for generation (default: {config.temperature})')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--mode', type=click.Choice(['api', 'local']), default=None,
              help=f'Execution mode: api (remote APIs) or local (default: {config.execution_mode})')
def cli(
    input_file: Path,
    output: Optional[Path],
    chunk_size: Optional[int],
    chunk_method: str,
    lines_per_chunk: int,
    provider: Optional[str],
    model: Optional[str],
    temperature: Optional[float],
    verbose: bool,
    mode: Optional[str]
) -> None:
    """
    Math Summarizer - chunk and summarize math textbooks using multiple LLM providers.
    
    Supports Groq, Together AI, Anthropic Claude, and Hugging Face models.
    
    INPUT_FILE: Path to the markdown file to process
    """
    if verbose:
        config.log_level = "DEBUG"
        logger.setLevel("DEBUG")
    
    # Override config with CLI options
    if chunk_size:
        config.max_chunk_size = chunk_size
    if provider:
        config.llm_provider = provider
        console.print(f"[blue]Using provider: {provider}[/blue]")
    if model:
        # Set the appropriate model based on provider
        if config.llm_provider == "groq":
            config.groq_model = model
        elif config.llm_provider == "together":
            config.together_model = model
        elif config.llm_provider == "anthropic":
            config.anthropic_model = model
        elif config.llm_provider == "huggingface":
            config.hf_model_name = model
        console.print(f"[blue]Using custom model: {model}[/blue]")
    if temperature is not None:
        config.temperature = temperature
    if mode is not None:
        config.execution_mode = mode
    
    # Set default output file
    if not output:
        output = input_file.parent / f"{input_file.stem}_summary.md"
    
    # Get current model name based on provider (using same logic as factory)
    from .llm_clients.factory import LLMClientFactory
    try:
        main_client = LLMClientFactory.create_main_client(config)
        current_model = main_client.model
        cleanup_model = config.latex_cleanup_model if hasattr(config, 'latex_cleanup_model') else "N/A"
    except Exception as e:
        # Fallback to old logic if factory fails
        if config.llm_provider == "groq":
            current_model = config.groq_model
        elif config.llm_provider == "together":
            current_model = config.together_model
        elif config.llm_provider == "anthropic":
            current_model = config.anthropic_model
        elif config.llm_provider == "huggingface":
            current_model = config.hf_model_name
        else:
            current_model = "unknown"
        cleanup_model = "N/A"
    
    console.print(f"[bold green]Math Summarizer v1.0.0[/bold green]")
    console.print(f"Input file: {input_file}")
    console.print(f"Output file: {output}")
    console.print(f"Provider: {config.llm_provider.upper()}")
    console.print(f"Main model: {current_model}")
    console.print(f"LaTeX cleanup model: {cleanup_model}")
    console.print(f"Execution mode: {config.execution_mode.upper()}")
    console.print(f"Chunk method: {chunk_method}")
    if chunk_method == 'tokens':
        console.print(f"Max chunk size: {config.max_chunk_size} tokens")
    else:
        console.print(f"Lines per chunk: {lines_per_chunk}")
    console.print()
    
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console,
        ) as progress:
            
            # Step 1: Chunk the document
            task1 = progress.add_task("Chunking document...", total=None)
            chunker = Chunker(config)
            
            if chunk_method == 'tokens':
                chunks = chunker.chunk_file(input_file)
            else:
                with open(input_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                chunks = chunker.chunk_by_lines(content, str(input_file), lines_per_chunk)
            
            progress.update(task1, description=f"Created {len(chunks)} chunks", completed=1, total=1)
            
            # Step 2: Summarize chunks
            task2 = progress.add_task(f"Summarizing chunks with {config.llm_provider.upper()} model...", total=len(chunks))
            summarizer = Summarizer(config)
            summaries = []
            
            for i, chunk in enumerate(chunks):
                summary = summarizer.summarize_chunk(chunk)
                summaries.append(summary)
                progress.update(task2, advance=1, 
                              description=f"Summarized {i+1}/{len(chunks)} chunks")
            
            # Step 3: Combine summaries
            task3 = progress.add_task("Combining summaries...", total=None)
            coalescer = Coalescer(config)
            final_summary = coalescer.coalesce_summaries(summaries, str(input_file))
            progress.update(task3, description="Combined summaries", completed=1, total=1)
        
        # Save output
        coalescer.save_summary(final_summary, output)
        
        console.print(f"\n[bold green]✓ Summary completed successfully![/bold green]")
        console.print(f"Output saved to: {output}")
        
        # Show some stats
        total_chars = sum(len(s) for s in summaries)
        console.print(f"\n[dim]Stats:[/dim]")
        console.print(f"[dim]  Chunks processed: {len(chunks)}[/dim]")
        console.print(f"[dim]  Total summary length: {total_chars:,} characters[/dim]")
        console.print(f"[dim]  Average per chunk: {total_chars // len(chunks):,} characters[/dim]")
        
    except Exception as e:
        console.print(f"\n[bold red]✗ Error: {str(e)}[/bold red]")
        logger.error(f"Processing failed: {str(e)}", exc_info=True)
        raise click.ClickException(str(e))


if __name__ == "__main__":
    cli()
