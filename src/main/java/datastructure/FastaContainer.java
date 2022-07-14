package datastructure;

/**
 * Represents the entry from a `.fasta` file.
 * <p>
 * Internal use to store the header and DNA sequence from a `.fasta` file.
 *
 * @author Alexander Seitz
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class FastaContainer {

    /**
     * Represents the header line from a `.fasta` entry.
     */
    private final String header;
    /**
     * Represents the DNA sequence from a `.fasta` entry.
     */
    private final char[] sequence;

    /**
     * Constructor of {@link FastaContainer}.
     *
     * @param header   {@link String} representing the header line of a `.fasta` entry.
     * @param sequence {@link Character[]} representing the DNA sequence of a `.fasta` entry.
     */
    public FastaContainer(String header, String sequence) {
        this.header = header;
        this.sequence = sequence.toCharArray();
    }

    /**
     * Returns the stored header line.
     *
     * @return {@link String} representing the header line of the represented `.fasta` entry.
     */
    public String getHeader() {
        return this.header;
    }

    /**
     * Returns the stored sequence in full length.
     *
     * @return {@link Character[]} representing the DNA sequence.
     */
    public String getSequence() {
        StringBuilder sequence = new StringBuilder();
        for (char c : this.sequence) {
            sequence.append(c);
        }
        return sequence.toString();
    }

    /**
     * Returns a sub sequence of the stored sequence.
     *
     * @param start The 1-based indexed starting position of the sequence to return.
     * @param end   The 1-based indexed end position of the sequence to return.
     * @return {@link Character[]} representing the DNA sequence from `start` to `end`.
     */
    public String getSequence(int start, int end) {
        if (start < 1 || end > this.sequence.length) {
            throw new IndexOutOfBoundsException("Accessed `.fasta` entry out of bounds of the stored sequence.");
        }
        StringBuilder sequence = new StringBuilder();
        for (int i = start; i <= end; i++) {
            sequence.append(this.sequence[i - 1]);
        }
        return sequence.toString();
    }
}