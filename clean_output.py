#!/usr/bin/env python3
import re

def clean_reasoning_traces(text):
    """Remove reasoning traces from the output."""
    original_length = len(text)
    
    # Remove <think> blocks completely
    text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL)
    
    # Remove other reasoning patterns
    reasoning_patterns = [
        r'Let me.*?(?=\n#|\n\n|$)',  # "Let me think about this..."
        r'I need to.*?(?=\n#|\n\n|$)',  # "I need to analyze..."
        r'First.*?(?=\n#|\n\n|$)',  # "First, let me..."
        r'Okay.*?(?=\n#|\n\n|$)',   # "Okay, so I need to..."
        r'Now.*?(?=\n#|\n\n|$)',    # "Now I'll process..."
        r'Starting with.*?(?=\n#|\n\n|$)',  # "Starting with section..."
        r'Similarly.*?(?=\n#|\n\n|$)',  # "Similarly, for section..."
    ]
    
    for pattern in reasoning_patterns:
        text = re.sub(pattern, '', text, flags=re.DOTALL)
    
    # Clean up multiple newlines
    text = re.sub(r'\n{3,}', '\n\n', text)
    text = text.strip()
    
    print(f"Cleaned: {original_length} -> {len(text)} characters ({100*(original_length-len(text))/original_length:.1f}% reduction)")
    return text

# Read the file
with open('output/ch1_summary.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Clean it
cleaned_content = clean_reasoning_traces(content)

# Save cleaned version
with open('output/ch1_summary_cleaned.md', 'w', encoding='utf-8') as f:
    f.write(cleaned_content)

print("Cleaned file saved as output/ch1_summary_cleaned.md")
