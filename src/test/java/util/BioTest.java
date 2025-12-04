package util;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.tuple.Triple;
import org.junit.jupiter.api.Test;
import uk.co.omegaprime.btreemap.BTreeMap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrowsExactly;

@SuppressWarnings("SpellCheckingInspection")
public class BioTest {

    @Test
    void globalNucleotideSequenceAlignmentIdenticalSequences() {
        String sequenceA = "TAAGTTTACA";
        String sequenceB = "TAAGTTTACA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAAGTTTACA", "TAAGTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentSingleSubstitution() {
        String sequenceA = "TAACTTTACA";
        String sequenceB = "TAAGTTTACA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAACTTTACA", "TAAGTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentSingleInsertion() {
        String sequenceA = "TAACTTTACA";
        String sequenceB = "TAATTTACA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAACTTTACA", "TAA-TTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentSingleDeletion() {
        String sequenceA = "TAATTTACA";
        String sequenceB = "TAACTTTACA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAA-TTTACA", "TAACTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentEmptySequenceA() {
        String sequenceA = "";
        String sequenceB = "TAATTTACA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("---------", "TAATTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentEmptySequenceB() {
        String sequenceA = "TAATTTACA";
        String sequenceB = "";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAATTTACA", "---------"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentSingelton() {
        String sequenceA = "A";
        String sequenceB = "G";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("A", "G"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentNoGapPrefix() {
        String sequenceA = "CC";
        String sequenceB = "CCCCCCC";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, false, 0);
        assertEquals(new Tuple<>("C-----C", "CCCCCCC"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentNoGapSuffix() {
        String sequenceA = "ATG";
        String sequenceB = "ATGCTACTTC";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, true, 0);
        assertEquals(new Tuple<>("A-------TG", "ATGCTACTTC"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentNoGapEnds() {
        String sequenceA = "CC";
        String sequenceB = "GCCCCG";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, true, 0);
        assertEquals(new Tuple<>("C----C", "GCCCCG"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentFreeGap() {
        String sequenceA = "AACA";
        String sequenceB = "TTATATCTA";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 0, 0,
                false, false, 0);
        assertEquals(new Tuple<>("--A-A-C-A", "TTATATCTA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentBanded() {
        String sequenceA = "TTTCGTATAACCTATGATAAAAAACTAACAATAATCATTAAATA";
        String sequenceB = "GCTGGATCGTATAACCAGCGGCGCCGCGCCTGGCCCACGGCTACCG";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 2);
        assertEquals(new Tuple<>(
                "TTTCGTATAACCTATGATAAAAAACTAACAATAATCA--TTAAATA",
                "GCTGGATCGTATAACCAGCGGCGCCGCGCCTGGCCCACGGCTACCG"
        ), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentComplexVariant() {
        String sequenceA = "CTGG";
        String sequenceB = "CCCCGAC";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, false, 0);
        assertEquals(new Tuple<>(
                "C---TGG",
                "CCCCGAC"
        ), result);
    }

    @Test
    void globalNucleotideSequenceAlignmentComplexFull() {
        String sequenceA = "CTAGACGCCGGGCCGCGGCCGTTGCCCATATTTAATATAAATTTTATCCCTACGGCGGCGCCGCATGCGGCCTCGGCGGC";
        String sequenceB = "TGCGTCACCCCCGCCCGCCCATATTTAATATAAATTTTATGCGACCCGCCCGAGAGGCGTGTATCGGGATCGGGTGGCGC";
        Tuple<String, String> result = Bio.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 8);
        assertEquals(new Tuple<>(
                "CTAGACGCCGGGCCGCGGCCGTTGCCCATATTTAATATAAATTTTAT----CCCTACGGCGGCGCCGCATGCGGCCTCG---GCGGC",
                "TGCGTCAC----CCCCGCCC---GCCCATATTTAATATAAATTTTATGCGACCCGCCCGAGAGGCGTGTATCGGGATCGGGTGGCGC"
        ), result);
    }

    @Test
    void globalProteinSequenceAlignmentIdenticalSequences() {
        String sequenceA = "MNLSVTLVRV";
        String sequenceB = "MNLSVTLVRV";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MNLSVTLVRV", "MNLSVTLVRV"), result);
    }

    @Test
    void globalProteinSequenceAlignmentSingleSubstitution() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVFGTIYSGLA", "MKETIPMQKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignmentInsertion() {
        String sequenceA = "MKETIPMQKNVEPAPYYFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVEPAPYYFGTIYSGLA", "MKETIPMQKNV------FGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignmentDeletion() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVEPAPYYFGTIYSGLA";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNV------FGTIYSGLA", "MKETIPMQKNVEPAPYYFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignmentIndel() {
        String sequenceA = "MKETIPMRTCEQQKNVFGTIYA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMRTCEQQKNVFGTIY---A", "MKETIPM-----QKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignmentEmptySequenceA() {
        String sequenceA = "";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("--------------------", "MKETIPMQKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignmentEmptySequenceB() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "";
        Tuple<String, String> result = Bio.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVFGTIYSGLA", "--------------------"), result);
    }

    @Test
    void getCanonicalVariantsIdenticalSequences() {
        String reference = "CGGGG";
        String alternative = "CGGGG";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(0, result.size());
    }

    @Test
    void getCanonicalVariantsSingleSubstitution() {
        String reference = "CGGGG";
        String alternative = "CGGTG";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(3, "G", "T"), result.get(0));
    }

    @Test
    void getCanonicalVariantsSingleInsertion() {
        String reference = "CGGGG-";
        String alternative = "CGGGGG";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(4, "G-", "GG"), result.get(0));
    }

    @Test
    void getCanonicalVariantsSingleDeletion() {
        String reference = "CGGGG";
        String alternative = "CGGG-";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(3, "GG", "G-"), result.get(0));
    }

    @Test
    void getCanonicalVariantsPolySubstitution() {
        String reference = "CGGGG";
        String alternative = "AGGGA";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(2, result.size());
        assertEquals(Triple.of(0, "C", "A"), result.get(0));
        assertEquals(Triple.of(4, "G", "A"), result.get(1));
    }

    @Test
    void getCanonicalVariantsComplex1() {
        String reference = "AGCTAGTCG---ATCTGCTAGT";
        String alternative = "AGCTAGTCGTTTATCTGCCAGT";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(2, result.size());
        assertEquals(Triple.of(8, "G---", "GTTT"), result.get(0));
        assertEquals(Triple.of(15, "T", "C"), result.get(1));
    }

    @Test
    void getCanonicalVariantsProteinSequences() {
        String reference =
                "MLKKASAFLIASCCVMSLAWAQANDNWYEGKPISAISFEGLEYIARGQLDTIFSQYKGQKWTYELYLEILQKVYDLEYFSEVSPKAVPTDPEYQYVMLQFTVKERPSVKGIKMVGNSQIRSGDLLSKILLKKGDIYNEVKMKVDQESLRRHYLDQGYAAVKISCEAKTEAGGVVVQFTIQEGKQTVVSRIQFKGNKAFTESVLKKVLSTQEARFLTSGVFKENALEADKAAVHSYYAERGYIDARVEGVAKTVDKKTDASRNLVTLTYTVVEGEQYRYGGVTIVGNQIFSTEELQAKIRLKRGAIMNMVAFEQGFQALADAYFENGYTSNYLNKEEHRDTAEKTLSFKITVVERERSHVEHIIIKGTKNTKDEVILREMLLKPGDVFSKSKFTDSLRNLFNLRYFSSLVPDVRPGSEQDLVDIILNVEEQSTANVQFGVTFSGVGEAGTFPLSLFCQWEEKNFLGKGNEISVNATLGSEAQSLKLGYVERWFLGSPLTVGFDFELTHKNLFVYRAGAKGNGLPHPYVSKEHWANSPGLAESFRLKYSRFESAIGAHTGYQWYPRYAVIRVNGGVDFRVVKNFYDKDNNQPFDLTVKEQLNWTSINSFWTSVSFDGRDFAYDPSSGWFLGQRCTFNGLVPCLEKEHSFRSDTKAEFYVTLLNYPVSAVWNLKFVLAFYTGVSVQTYYGRRKSENGKGNGVRSGALVIDGVLVGRGWSEDAKKNTGDLLLHHWIEFRWPLAHGIVSFDFFFDAAMVYNIESQSPNGSSSASSSSSSSSSSSSTTSS----EGLYKMSYGPGLRFTLPQFPLKLAFANTFTSPGGIPKTKKDWNFVLSFTVNNL";
        String alternative =
                "MLKKASAFLIASCCVMSLAWAQANDNWYEGKPISAISFEGLEYIARGQLDTIFSQYKGQKWTYELYLEILQKVYDLEYFSEVSPKAVPTDPEYQYVMLQFTVKERPSVKGIKMVGNSQIRSGDLLSKILLKKGDIYNEVKMKVDQESLRRHYLDQGYAAVKISCEAKTEAGGVVVQFTIQEGKQTVVSRIQFKGNKAFTESVLKKVLSTQEARFLTSGVFKENALEADKAAVHSYYAERGYIDARVEGVAKTVDKKTDASRNLVTLTYTVVEGEQYRYGGVTIVGNQIFSTEELQAKIRLKRGAIMNMVAFEQGFQALADAYFENGYTSNYLNKEEHRDTAEKTLSFKITVVERERSHVEHIIIKGTKNTKDEVILREMLLKPGDVFSKSKFTDSLRNLFNLRYFSSLVPDVRPGSEQDLVDIILNVEEQSTANVQFGVTFSGVGEAGTFPLSLFCQWEEKNFLGKGNEISVNATLGSEAQSLKLGYVERWFLGSPLTVGFDFELTHKNLFVYRAGSYGNGLPHPYTSREQWASSPGLAESFRLKYSRFESAIGAHTGYQWYPRYAVIRVNGGVDFRVVKNFYDKDNNQPFDLTVEEQLNWTSINSFWTSVSFDGRDFAYDPSSGWFLGQRCTFNGLVPFLEKEHSFRSDTKAEFYVTLLNYPVSAVWNLKFVLAFYTGVSVQTYYGRRKSENGKGNGVRSGALVIDGVLVGRGWSEDAKKNTGDLLLHHWIEFRWPLAHGIVSFDFFFDAAMVYNIESQSPNGSSSASSSSSSSSSSSSSSSSSSSSEGLYKMSYGPGLRFTLPQFPLKLAFANTFTSPGGIPKTKKNWNFVLSFTVNNL";
        ArrayList<Triple<Integer, String, String>> result = Bio.getCanonicalVariants(reference, alternative);
        assertEquals(12, result.size());
        assertEquals(Triple.of(516, "A", "S"), result.get(0));
        assertEquals(Triple.of(517, "K", "Y"), result.get(1));
        assertEquals(Triple.of(526, "V", "T"), result.get(2));
        assertEquals(Triple.of(528, "K", "R"), result.get(3));
        assertEquals(Triple.of(530, "H", "Q"), result.get(4));
        assertEquals(Triple.of(533, "N", "S"), result.get(5));
        assertEquals(Triple.of(595, "K", "E"), result.get(6));
        assertEquals(Triple.of(639, "C", "F"), result.get(7));
        assertEquals(Triple.of(780, "T", "S"), result.get(8));
        assertEquals(Triple.of(781, "T", "S"), result.get(9));
        assertEquals(Triple.of(783, "S----", "SSSSS"), result.get(10));
        assertEquals(Triple.of(824, "D", "N"), result.get(11));
    }

    @Test
    void getCanonicalVariantsException() {
        String reference = "GTATGGGGCT";
        String alternative = "GGGGGCT";
        assertThrowsExactly(IllegalArgumentException.class, () -> Bio.getCanonicalVariants(reference, alternative));
    }

    @Test
    void translateSequenceValidSequence() throws MusialException {
        String sequence = "ATGCGT";
        String result = Bio.translateSequence(sequence, false);
        assertEquals("MR", result);
    }

    @Test
    void translateSequenceReverseSequence() throws MusialException {
        String sequence = "ATGCGT";
        String result = Bio.translateSequence(sequence, true);
        assertEquals("TH", result);
    }

    @Test
    void translateSequenceEmptySequence() throws MusialException {
        String sequence = "";
        String result = Bio.translateSequence(sequence, false);
        assertEquals("", result);
    }

    @Test
    void translateSequenceInvalidSequence() {
        String sequence = "ATGCGTX";
        assertThrowsExactly(MusialException.class, () -> Bio.translateSequence(sequence, false));
    }

    @Test
    void integrateVariantsByStringEmptyReference() {
        String reference = "";
        Map<Integer, String> variants = new HashMap<>();
        variants.put(0, "A");
        assertThrowsExactly(IllegalArgumentException.class, () -> Bio.integrateVariants(reference, variants));
    }

    @Test
    void integrateVariantsByStringEmptyVariants() throws MusialException {
        String reference = "ACGTACGT";
        Map<Integer, String> variants = new HashMap<>();
        String result = Bio.integrateVariants(reference, variants);
        assertEquals(reference, result);
    }

    @Test
    void integrateVariantsByStringNonCanonicalDeletion() {
        String reference = "ACGTACGT";
        Map<Integer, String> variants = new HashMap<>();
        variants.put(2, "---");
        assertThrowsExactly(MusialException.class, () -> Bio.integrateVariants(reference, variants));
    }

    @Test
    void integrateVariantsByStringNonCanonicalInDel() {
        String reference = "ACGTACGT";
        Map<Integer, String> variants = new HashMap<>();
        variants.put(2, "A--TG");
        assertThrowsExactly(MusialException.class, () -> Bio.integrateVariants(reference, variants));
    }

    @Test
    void integratesVariantsByString1() throws MusialException {
        String reference = "ACGTACGT";
        Map<Integer, String> variants = new HashMap<>();
        variants.put(2, "T");
        variants.put(5, "G-");
        String result = Bio.integrateVariants(reference, variants);
        assertEquals("ACTTAGT", result);
    }

    @Test
    void integratesVariantsByString2() throws MusialException {
        String reference = "ACGTACGT";
        Map<Integer, String> variants = new HashMap<>();
        variants.put(2, "T-");
        variants.put(3, "A");
        variants.put(7, "TGAT");
        String result = Bio.integrateVariants(reference, variants);
        assertEquals("ACTACGTGAT", result);
    }

    @Test
    void integrateVariantsEmptyReference() {
        Bio.ReferenceContext[] reference = new Bio.ReferenceContext[0];
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(1, "T");
        assertThrowsExactly(IllegalArgumentException.class, () -> Bio.integrateVariants(reference, variants, false));
    }

    @Test
    void integrateVariantsEmptyVariants() throws MusialException {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 0),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(6, 'C', 0),
                new Bio.ReferenceContext(7, 'G', 0),
                new Bio.ReferenceContext(8, 'T', 0)
        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        String result = Bio.integrateVariants(reference, variants, false);
        assertEquals("ACGTACGT", result);
    }

    @Test
    void integrateVariantsNonCanonicalDeletion() {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 0),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(6, 'C', 0),
                new Bio.ReferenceContext(7, 'G', 0),
                new Bio.ReferenceContext(8, 'T', 0)
        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(2, "---");
        assertThrowsExactly(MusialException.class, () -> Bio.integrateVariants(reference, variants, false));
    }

    @Test
    void integrateVariantsNonCanonicalInDel() {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 0),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(6, 'C', 0),
                new Bio.ReferenceContext(7, 'G', 0),
                new Bio.ReferenceContext(8, 'T', 0)
        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(2, "C--TG");
        assertThrowsExactly(MusialException.class, () -> Bio.integrateVariants(reference, variants, false));
    }

    @Test
    void integrateVariantsFaultyReference() {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 0),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(8, 'T', 0),
                new Bio.ReferenceContext(6, 'C', 0),
                new Bio.ReferenceContext(7, 'G', 0)
        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(1, "G");
        assertThrowsExactly(MusialException.class, () -> Bio.integrateVariants(reference, variants, false));
    }

    @Test
    void integratesVariants1() throws MusialException {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 1),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(6, 'C', 1),
                new Bio.ReferenceContext(7, 'G', 4),
                new Bio.ReferenceContext(8, 'T', 1)

        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(1, "G");
        variants.put(4, "T--");
        variants.put(7, "GAT");
        String result = Bio.integrateVariants(reference, variants, false);
        assertEquals("GCGT----GAT--T-", result);
    }

    @Test
    void integratesVariants2() throws MusialException {
        Bio.ReferenceContext[] reference = {
                new Bio.ReferenceContext(1, 'A', 0),
                new Bio.ReferenceContext(2, 'C', 0),
                new Bio.ReferenceContext(3, 'G', 0),
                new Bio.ReferenceContext(4, 'T', 1),
                new Bio.ReferenceContext(5, 'A', 0),
                new Bio.ReferenceContext(6, 'C', 1),
                new Bio.ReferenceContext(7, 'G', 5),
                new Bio.ReferenceContext(8, 'T', 1)

        };
        NavigableMap<Integer, String> variants = BTreeMap.create();
        variants.put(1, "G");
        variants.put(4, "T--");
        variants.put(7, "GAT");
        String result = Bio.integrateVariants(reference, variants, true);
        assertEquals("GCGTGATT", result);
    }

    @Test
    void alignByCigarMatchOnly() {
        String reference = "ACGTACGT";
        String query = "ACGTACGT";
        String cigar = "8M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGTACGT", "ACGTACGT"), result);
    }

    @Test
    void alignByCigarInsertion() {
        String reference = "ACGTACGT";
        String query = "ACGTTACGT";
        String cigar = "4M1I4M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGT-ACGT", "ACGTTACGT"), result);
    }

    @Test
    void alignByCigarDeletion() {
        String reference = "ACGTACGT";
        String query = "ACGACGT";
        String cigar = "3M1D4M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGTACGT", "ACG-ACGT"), result);
    }

    @Test
    void alignByCigarSoftClipping() {
        String reference = "ACGTACGT";
        String query = "TTACGTACGT";
        String cigar = "2S8M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGTACGT", "ACGTACGT"), result);
    }

    @Test
    void alignByCigarHardClipping() {
        String reference = "ACGTACGT";
        String query = "ACGTACGT";
        String cigar = "2H8M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGTACGT", "ACGTACGT"), result);
    }

    @Test
    void alignByCigarComplexCigar() {
        String reference = "ACGTACGTA";
        String query = "ACGTTACGA";
        String cigar = "4M1I3M1D1M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 0);
        assertEquals(new Tuple<>("ACGT-ACGTA", "ACGTTACG-A"), result);
    }

    @Test
    void alignByCigarUnsupportedOperation() {
        String reference = "ACGTACGT";
        String query = "ACGTACGT";
        String cigar = "8Z";
        assertThrowsExactly(IllegalArgumentException.class, () -> Bio.alignByCigar(reference, query, cigar, 0));
    }

    @Test
    void alignByCigarMismatchedLengths() {
        String reference = "ACGTACGT";
        String query = "ACGT";
        String cigar = "8M";
        assertThrowsExactly(IllegalArgumentException.class, () -> Bio.alignByCigar(reference, query, cigar, 0));
    }

    @Test
    void alignByCigarOffsetHandling() {
        String reference = "ACGTACGT";
        String query = "GTACGT";
        String cigar = "6M";
        Tuple<String, String> result = Bio.alignByCigar(reference, query, cigar, 2);
        assertEquals(new Tuple<>("GTACGT", "GTACGT"), result);
    }

}
