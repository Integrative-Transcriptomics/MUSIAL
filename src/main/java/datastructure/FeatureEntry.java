package datastructure;

import exceptions.MusialBioException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.javatuples.Triplet;
import components.Bio;
import components.IO;
import components.Logging;

/**
 * Internal representation of a reference sequence location that is subject to analysis. May represent the full
 * genome, a single gene, contigs or plasmids and chromosomes.
 * <p>
 * Pairs of instances of this class and {@link SampleEntry} instances are used internally as so called 'Run
 * entries' to specify a set of analysis tasks, i.e. which sample is analyzed with respect to which specified
 * reference feature.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.0
 */
public final class FeatureEntry {

    /**
     * The internal name of the entry.
     */
    public final String name;
    /**
     * The reference DNA sequence of the entry.
     */
    public String nucleotideSequence;
    /**
     * The location of the entry on the reference, i.e. for genes the contig or chromosome the feature is located on.
     */
    public final String chromosome;
    /**
     * The 1-based indexed starting position of the feature.
     */
    public final int start;
    /**
     * The 1-based indexed end position of the feature.
     */
    public final int end;
    /**
     * Indicates if the feature is located on the sense strand.
     */
    public final boolean isSense;
    /**
     * A {@link AllocatedProteinEntry} object containing information about the reference protein structure of the possible
     * protein product of the gene represented by this {@link FeatureEntry} and all discovered proteoforms.
     */
    public AllocatedProteinEntry allocatedProtein;
    /**
     * A {@link File} object pointing to a .pdb format file containing the protein structure information of
     * {@link FeatureEntry#allocatedProtein}.
     */
    public transient File pdbFile;
    /**
     * A {@link HashMap} yielding any meta-information {@link String} key-value pairs about this {@link FeatureEntry}.
     */
    public final HashMap<String, String> annotations = new HashMap<>();

    /**
     * Constructor of {@link FeatureEntry}.
     *
     * @param entryName     {@link String} representing the internal name of the reference feature to analyze.
     * @param entryLocation {@link String} the name of the reference location (contig, chromosome, plasmid) the feature
     *                      is located on.
     * @param entryStart    {@link Integer} The 1-based indexed starting position of the feature on the reference.
     * @param entryEnd      {@link Integer} The 1-based indexed end position of the feature on the reference.
     * @throws MusialBioException If the specified locus is ambiguous.
     */
    public FeatureEntry(String entryName, String entryLocation,
                        int entryStart,
                        int entryEnd) throws MusialBioException {
        this.name = entryName;
        this.chromosome = entryLocation;
        if (entryStart >= 0 && entryEnd > 0 && (entryEnd > entryStart)) {
            // CASE: Feature is on sense strand.
            this.isSense = true;
            this.start = entryStart + 1;
            this.end = entryEnd;
        } else if (entryStart < 0 && entryEnd <= 0 && (entryEnd < entryStart)) {
            // CASE: Feature is on anti-sense strand.
            this.isSense = false;
            this.start = -entryStart + 1;
            this.end = -entryEnd;
        } else {
            throw new MusialBioException(
                    "Failed to add gene feature due to faulty position data (start, end)\t" + entryStart + ", " + entryEnd);
        }
    }

    /**
     * Method to add the nucleotide sequence to any {@link FeatureEntry}.
     *
     * @param featureEntry      The {@link FeatureEntry} to which nucleotide sequence information shall be added.
     * @param referenceSequence {@link FastaContainer} representation of the parent/reference sequence, i.e. reference
     *                          genome or genome assembly.
     */
    public static void imputeSequence(FeatureEntry featureEntry, FastaContainer referenceSequence) {
        featureEntry.nucleotideSequence = referenceSequence.getSequence(featureEntry.start, featureEntry.end);
    }

    /**
     * Method to add the protein wild-type information to any {@link FeatureEntry}.
     * <p>
     * The specified {@link FeatureEntry} has to specify the {@link FeatureEntry#pdbFile} property.
     * <p>
     * - The parsed protein amino-acid sequence is compared with the translated amino-acid sequence of the
     * {@link FeatureEntry#nucleotideSequence} property by pairwise sequence alignment. While the translated sequence
     * is allowed to be longer than the protein sequence, this is not allowed vice versa.
     * - The internally stored {@link FeatureEntry#allocatedProtein}s {@link AllocatedProteinEntry#pdb} will be matched
     * in terms of residue numbers to the information gained from the alignment, i.e. if the .pdb starts with residue 10,
     * but the sequences align at the beginning the first residue number in the internal .pdb is set to be 1; if the
     * .pdb starts with residue 1, but the sequence (of the protein) has a gap of length 9 at the beginning in the
     * alignment the first residue number in the internal .pdb is set to be 10.
     *
     * @param featureEntry The {@link FeatureEntry} to which protein information shall be added.
     * @throws IOException        If the .pdb file pointed to by the {@link FeatureEntry#pdbFile} can not be parsed.
     * @throws MusialBioException If either the translated {@link FeatureEntry#nucleotideSequence} yields internal
     *                            termination codons or contains any gaps when aligned to the amino-acid sequence derived from the .pdb fiel.
     */
    public static void imputeProtein(FeatureEntry featureEntry) throws IOException, MusialBioException {
        Structure pdbStructure = IO.readStructure(featureEntry.pdbFile);
        HashMap<String, String> proteinSequences = IO.getSequencesFromPdbStructure(pdbStructure);
        String chainId;
        String chainSeq;
        char[] alignedChainSeq;
        String translatedFeatureSequence = Bio.translateNucSequence(featureEntry.nucleotideSequence, true, true,
                featureEntry.isSense);
        if (translatedFeatureSequence.endsWith(String.valueOf(Bio.TERMINATION_AA1))) {
            if (translatedFeatureSequence.substring(0, translatedFeatureSequence.length() - 1).contains(String.valueOf(Bio.TERMINATION_AA1))) {
                throw new MusialBioException("Failed to allocate protein " + featureEntry.pdbFile.getName() + " to " +
                        "feature " + featureEntry.name + " due to internal translated termination in translated feature nucleotide sequence.");
            }
        } else {
            Logging.logWarning("Feature " + featureEntry.name + " does not end with a translated termination and may be inappropriate.");
        }
        char[] alignedTranslFeatureSequence;
        StringBuilder paddedChainSeqBuilder = new StringBuilder();
        List<String> diffSegments;
        for (Map.Entry<String, String> proteinSequenceEntry : proteinSequences.entrySet()) {
            paddedChainSeqBuilder.setLength(0);
            chainId = proteinSequenceEntry.getKey();
            chainSeq = proteinSequenceEntry.getValue();
            // Align protein and transl. reference sequence to infer correct starting position and truncated prefix peptides.
            Triplet<Integer, String, String> alignedSequences =
                    Bio.globalAminoAcidSequenceAlignment(
                            chainSeq,
                            translatedFeatureSequence,
                            5,
                            4,
                            Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FREE,
                            Bio.GLOBAL_SEQUENCE_ALIGNMENT_MARGIN_GAP_MODES.FREE
                    );
            alignedChainSeq = alignedSequences.getValue1().toCharArray();
            alignedTranslFeatureSequence = alignedSequences.getValue2().toCharArray();
            for (int pos = 0; pos < alignedChainSeq.length; pos++) {
                if (alignedChainSeq[pos] != alignedTranslFeatureSequence[pos]) {
                    if (alignedChainSeq[pos] == Bio.GAP) {
                        paddedChainSeqBuilder.append(Character.toLowerCase(alignedTranslFeatureSequence[pos]));
                    } else if (alignedTranslFeatureSequence[pos] == Bio.GAP) {
                        throw new MusialBioException("Failed to allocate protein " + featureEntry.pdbFile.getName() + " to " +
                                "feature " + featureEntry.name + " due to gaps in the aligned translated feature nucleotide sequence.");
                    }
                } else {
                    paddedChainSeqBuilder.append(alignedChainSeq[pos]);
                }
            }
            // Construct Iterator of current chain groups sequence numbers.
            Chain pdbChain = pdbStructure.getChain(chainId);
            List<Group> pdbChainAtomGroups = pdbChain.getAtomGroups();
            Iterator<Group> pdbChainAtomGroupIterator = pdbChainAtomGroups.iterator();
            List<Group> pdbChainAtomGroupsFixed = new ArrayList<>();
            char[] paddedChainSeqChars = paddedChainSeqBuilder.toString().toCharArray();
            for (int chainPosition = 1; chainPosition < paddedChainSeqChars.length + 1; chainPosition++) {
                if (Character.isLowerCase(paddedChainSeqChars[chainPosition - 1])) {
                    continue;
                }
                if (pdbChainAtomGroupIterator.hasNext()) {
                    Group nextGroup = pdbChainAtomGroupIterator.next();
                    nextGroup.getResidueNumber().setSeqNum(chainPosition);
                    pdbChainAtomGroupsFixed.add(nextGroup);
                } else {
                    break;
                }
            }
            pdbChain.setAtomGroups(pdbChainAtomGroupsFixed);
            diffSegments = Arrays.stream(
                    paddedChainSeqBuilder.toString().split("(?=\\p{Upper})")
            ).filter(s -> s.length() > 1).collect(Collectors.toList());
            if (diffSegments.size() > 2) {
                Logging.logWarning("Feature " + featureEntry.name + " disaccords in " + diffSegments.size() + " segments " +
                        "with allocated protein " + featureEntry.pdbFile.getName() + " chain " + chainId + ": The structure may be inappropriate.");
            }
            proteinSequences.put(chainId, paddedChainSeqBuilder.toString());
        }
        featureEntry.allocatedProtein = new AllocatedProteinEntry(
                FilenameUtils.removeExtension(featureEntry.pdbFile.getName()),
                String.join(IO.LINE_SEPARATOR, pdbStructure.toPDB()),
                proteinSequences
        );
    }

}