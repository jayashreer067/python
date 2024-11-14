import random
import time

from Bio.Seq import Seq
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor

from concurrent.futures import ProcessPoolExecutor, TimeoutError
import logging 


def longest_contiguous(seq):
    """
    Count the length of the longest contiguous subsequences.

    This function finds each contiguous subsequence in ``seq``,
    e.g. ``AAAA`` or ``BB`` and counts its length and returns the
    letter with the longest contiguous subsequence and its length.

    It returns a dictionary with the key being the letter and the
    value being the length of the longest subsequence. If two letters
    have longest subsequences of the same length, both are included
    separately in the dictionary.

    Each unique letter in ``seq`` is counted separately so the
    sequence "AABBCC" will find the longest subsequences of each of
    "A", "B" and "C".

    If ``seq`` is empty then an empty dictionary is returned.

    Args:
        seq (list of Seq): A list of sequences.

    Returns:
        dict: The count of the longest subsequences, keyed by letter.
        Examples:
        >>> longest_contiguous("aacbbb")
        {'b': 3}
        >>> longest_contiguous("aabbbaabbc")
        {'b': 3}
        >>> longest_contiguous("aaabbbaabb")
        {'a': 3, 'b': 3}
        >>> longest_contiguous("")
        {}
    """
    # Check if seq is empty
    if not seq:
        return {}
     # Check if seq is an empty string
    if isinstance(seq, str) and not seq:
        raise ValueError("Input 'seq' cannot be an empty string")

    # Initialize an empty dictionary to store counts of longest subsequences
    longest_subsequences = {}

    # Check if seq is a string or Bio.Seq object
    if isinstance(seq[0], str):
        max_lengths = {}
        current_length = 1
        # for looping over each letter
        for i in range(1, len(seq)):
            if seq[i] == seq[i - 1]:
                current_length += 1
            else:
                if seq[i - 1] not in max_lengths or current_length > max_lengths[seq[i - 1]]:
                    max_lengths[seq[i - 1]] = current_length
                current_length = 1

        # Check the last character
        if seq[-1] not in max_lengths or current_length > max_lengths[seq[-1]]:
            max_lengths[seq[-1]] = current_length
            
        # Filter the dictionary to keep only the characters with the maximum subsequence length     
        longest_subsequences = {letter: length for letter, length in max_lengths.items() if length == max(max_lengths.values())}
        return longest_subsequences
    #if the seq is Bio.Seq object
    else:
        longest_subsequences = {}
        # for looping over each Seq object
        for i, sequence in enumerate(seq):
            seq_string = sequence.__str__()
            longest_subsequence = {}
            current_length = 0
            prev = None
            
            # for looping over each letter in the sequence
            for char in seq_string:
                if char == prev:
                    current_length += 1
                else:
                    current_length = 1
                    prev = char
                    
                # Check if the current character is already in the longest subsequence dictionary
                if char in longest_subsequence:
                    longest_subsequence[char] = max(longest_subsequence[char], current_length)
                else:
                    longest_subsequence[char] = current_length
                    
            #to get the maximum subsequence length among all the values   
            max_length = max(longest_subsequence.values())
            # Filter the dictionary to keep only the characters with the maximum subsequence length
            longest_subsequence = {char: length for char, length in longest_subsequence.items() if length == max_length}
            # Add the longest subsequence for the current sequence to the main dictionary
            longest_subsequences[i] = longest_subsequence
        # Return the dictionary containing longest subsequences for each sequence
        return longest_subsequences


    ...



def all_longest(each_longest):
    """
    Find the longest subsequences across multiple sequences.

    This function combines the outputs from ``longest_contiguous``
    to a single dictionary. For each letter in each of the input
    dictionaries, the output will contain that letter as a key and
    the largest of all associated values as its value.

    Args:
        each_longest (list): A list of longest subsequence counts.

    Returns:
        dict: The largest values for each letter from the inputs.

    Examples:
        >>> all_longest([{"a": 4}, {"b": 6}])
        {'a': 4, 'b': 6}
        >>> all_longest([{"a": 4}, {"b": 6, "a": 2}])
        {'a': 4, 'b': 6}
        >>> all_longest([{"a": 4}, {"b": 10, "a": 2}])
        {'a': 4, 'b': 10}
        >>> all_longest([{"a": 4}, {"b": 6, "a": 2}, {"b": 10}])
        {'a': 4, 'b': 10}
    """
    # Initialize an empty dictionary to store the longest subsequences
    combined_longest = {}

    # Iterate over each dictionary in each_longest
    for subsequence in each_longest:
        # Iterate over each letter and its length in the current subsequence
        for letter, length in subsequence.items():
                        # Check if length is numeric
            if not isinstance(length, (int, float)):
                raise ValueError("Values in dictionaries must be numeric.")
            # Update combined_longest with the maximum length for each letter
            combined_longest[letter] = max(combined_longest.get(letter, 0), length)

    return combined_longest

    ...


def generate_random_test_data(num_seq=100_000, seq_length=1000):
    """
    Generate a large amount of example data to test the above functions.

    Args:
        num_seq (int): The number of sequnces to return.
        seq_length (int): The length of each individual sequence.
    """
    sequences = []
    for i in range(num_seq):
        sequences.append(Seq("".join(random.choices(["g", "a", "t", "c"], k=seq_length ))))
    return sequences


def run_large_test_serial(sequences):
    """
    Run through a list of provided sequences in series and return the result.
    """
    longest_each_sequence = []
    for seq in sequences:
        longest_each_sequence.append(longest_contiguous(seq))

    answer = all_longest(longest_each_sequence)

    return answer



def run_large_test_parallel(sequences):
    """
    Run through a list of provided sequences in parallel and return the result.
    
    Args:
    - sequences: A list of sequences to be processed
    
    Returns:
    - The combined result of processing all sequences
    """
    try:
        # Create a ProcessPoolExecutor with a maximum of 10 worker processes
        with ProcessPoolExecutor(max_workers=10) as pool:
            # Submit each sequence to the executor for processing
            # Using the map method to process sequences in parallel
            # chunksize parameter specifies the number of items to send to each worker process at a time
            results = pool.map(longest_contiguous, sequences, chunksize=1000)

        # Combine the results obtained from parallel processing
        answer = all_longest(results)

        return answer

    except TimeoutError as e:
        logging.error(f"TimeoutError occurred: {e}")
        return None  # Return None to indicate failure

    except Exception as e:
        logging.error(f"An error occurred during parallel processing: {e}")
        return None  # Return None to indicate failure

if __name__ == "__main__":
    print("Generating test sequences")
    sequences = generate_random_test_data()

    print("Running serial test")
    start = time.perf_counter()  # start the performance timer
    answer1 = run_large_test_serial(sequences)  # Run the analysis serially
    duration_seconds = time.perf_counter() - start  # stop the performance timer
    print(f"  Serial: {duration_seconds:.2f} seconds - {answer1}")

    print("Running parallel test")
    start = time.perf_counter()
    answer2 = run_large_test_parallel(sequences)  # Run the parallel analysis
    duration_seconds = time.perf_counter() - start
    print(f"Parallel: {duration_seconds:.2f} seconds - {answer2}")
