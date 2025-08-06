package utility;

import exceptions.MusialException;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.tuple.Triple;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrowsExactly;

@SuppressWarnings("SpellCheckingInspection")
public class SequenceOperationsTest {

    @Test
    void globalNucleotideSequenceAlignment_identicalSequences() {
        String sequenceA = "TAAGTTTACA";
        String sequenceB = "TAAGTTTACA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAAGTTTACA", "TAAGTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleSubstitution() {
        String sequenceA = "TAACTTTACA";
        String sequenceB = "TAAGTTTACA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAACTTTACA", "TAAGTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleInsertion() {
        String sequenceA = "TAACTTTACA";
        String sequenceB = "TAATTTACA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAACTTTACA", "TAA-TTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleDeletion() {
        String sequenceA = "TAATTTACA";
        String sequenceB = "TAACTTTACA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAA-TTTACA", "TAACTTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_emptySequenceA() {
        String sequenceA = "";
        String sequenceB = "TAATTTACA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("---------", "TAATTTACA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_emptySequenceB() {
        String sequenceA = "TAATTTACA";
        String sequenceB = "";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("TAATTTACA", "---------"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singelton() {
        String sequenceA = "A";
        String sequenceB = "G";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 0);
        assertEquals(new Tuple<>("A", "G"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_noGapPrefix() {
        String sequenceA = "CC";
        String sequenceB = "CCCCCCC";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, false, 0);
        assertEquals(new Tuple<>("C-----C", "CCCCCCC"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_noGapSuffix() {
        String sequenceA = "ATG";
        String sequenceB = "ATGCTACTTC";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, true, 0);
        assertEquals(new Tuple<>("A-------TG", "ATGCTACTTC"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_noGapEnds() {
        String sequenceA = "CC";
        String sequenceB = "GCCCCG";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, true, 0);
        assertEquals(new Tuple<>("C----C", "GCCCCG"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_freeGap() {
        String sequenceA = "AACA";
        String sequenceB = "TTATATCTA";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 0, 0,
                false, false, 0);
        assertEquals(new Tuple<>("--A-A-C-A", "TTATATCTA"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_banded() {
        String sequenceA = "TTTCGTATAACCTATGATAAAAAACTAACAATAATCATTAAATA";
        String sequenceB = "GCTGGATCGTATAACCAGCGGCGCCGCGCCTGGCCCACGGCTACCG";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 2);
        assertEquals(new Tuple<>(
                "TTTCGTATAACCTATGATAAAAAACTAACAATAATCA--TTAAATA",
                "GCTGGATCGTATAACCAGCGGCGCCGCGCCTGGCCCACGGCTACCG"
        ), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_complexVariant() {
        String sequenceA = "CTGG";
        String sequenceB = "CCCCGAC";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                true, false, 0);
        assertEquals(new Tuple<>(
                "C---TGG",
                "CCCCGAC"
        ), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_complexFull() {
        String sequenceA = "CTAGACGCCGGGCCGCGGCCGTTGCCCATATTTAATATAAATTTTATCCCTACGGCGGCGCCGCATGCGGCCTCGGCGGC";
        String sequenceB = "TGCGTCACCCCCGCCCGCCCATATTTAATATAAATTTTATGCGACCCGCCCGAGAGGCGTGTATCGGGATCGGGTGGCGC";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 4, 2,
                false, false, 8);
        assertEquals(new Tuple<>(
                "CTAGACGCCGGGCCGCGGCCGTTGCCCATATTTAATATAAATTTTAT----CCCTACGGCGGCGCCGCATGCGGCCTCG---GCGGC",
                "TGCGTCAC----CCCCGCCC---GCCCATATTTAATATAAATTTTATGCGACCCGCCCGAGAGGCGTGTATCGGGATCGGGTGGCGC"
        ), result);
    }

    @Test
    void globalProteinSequenceAlignment_identicalSequences() {
        String sequenceA = "MNLSVTLVRV";
        String sequenceB = "MNLSVTLVRV";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MNLSVTLVRV", "MNLSVTLVRV"), result);
    }

    @Test
    void globalProteinSequenceAlignment_singleSubstitution() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVFGTIYSGLA", "MKETIPMQKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignment_insertion() {
        String sequenceA = "MKETIPMQKNVEPAPYYFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVEPAPYYFGTIYSGLA", "MKETIPMQKNV------FGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignment_deletion() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "MKETIPMQKNVEPAPYYFGTIYSGLA";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNV------FGTIYSGLA", "MKETIPMQKNVEPAPYYFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignment_indel() {
        String sequenceA = "MKETIPMRTCEQQKNVFGTIYA";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMRTCEQQKNVFGTIY---A", "MKETIPM-----QKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignment_emptySequenceA() {
        String sequenceA = "";
        String sequenceB = "MKETIPMQKNVFGTIYSGLA";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("--------------------", "MKETIPMQKNVFGTIYSGLA"), result);
    }

    @Test
    void globalProteinSequenceAlignment_emptySequenceB() {
        String sequenceA = "MKETIPMQKNVFGTIYSGLA";
        String sequenceB = "";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 15, 3,
                false, false, 0);
        assertEquals(new Tuple<>("MKETIPMQKNVFGTIYSGLA", "--------------------"), result);
    }

    @Test
    void getCanonicalVariants_identicalSequences() {
        String reference = "CGGGG";
        String alternative = "CGGGG";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(0, result.size());
    }

    @Test
    void getCanonicalVariants_singleSubstitution() {
        String reference = "CGGGG";
        String alternative = "CGGTG";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(3, "G", "T"), result.get(0));
    }

    @Test
    void getCanonicalVariants_singleInsertion() {
        String reference = "CGGGG-";
        String alternative = "CGGGGG";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(4, "G-", "GG"), result.get(0));
    }

    @Test
    void getCanonicalVariants_singleDeletion() {
        String reference = "CGGGG";
        String alternative = "CGGG-";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(1, result.size());
        assertEquals(Triple.of(3, "GG", "G-"), result.get(0));
    }

    @Test
    void getCanonicalVariants_polySubstitution() {
        String reference = "CGGGG";
        String alternative = "AGGGA";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(2, result.size());
        assertEquals(Triple.of(0, "C", "A"), result.get(0));
        assertEquals(Triple.of(4, "G", "A"), result.get(1));
    }

    @Test
    void getCanonicalVariants_complex1() {
        String reference = "AGCTAGTCG---ATCTGCTAGT";
        String alternative = "AGCTAGTCGTTTATCTGCCAGT";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
        assertEquals(2, result.size());
        assertEquals(Triple.of(8, "G---", "GTTT"), result.get(0));
        assertEquals(Triple.of(15, "T", "C"), result.get(1));
    }

    @Test
    void getCanonicalVariants_proteinSequences() {
        String reference = "MLKKASAFLIASCCVMSLAWAQANDNWYEGKPISAISFEGLEYIARGQLDTIFSQYKGQKWTYELYLEILQKVYDLEYFSEVSPKAVPTDPEYQYVMLQFTVKERPSVKGIKMVGNSQIRSGDLLSKILLKKGDIYNEVKMKVDQESLRRHYLDQGYAAVKISCEAKTEAGGVVVQFTIQEGKQTVVSRIQFKGNKAFTESVLKKVLSTQEARFLTSGVFKENALEADKAAVHSYYAERGYIDARVEGVAKTVDKKTDASRNLVTLTYTVVEGEQYRYGGVTIVGNQIFSTEELQAKIRLKRGAIMNMVAFEQGFQALADAYFENGYTSNYLNKEEHRDTAEKTLSFKITVVERERSHVEHIIIKGTKNTKDEVILREMLLKPGDVFSKSKFTDSLRNLFNLRYFSSLVPDVRPGSEQDLVDIILNVEEQSTANVQFGVTFSGVGEAGTFPLSLFCQWEEKNFLGKGNEISVNATLGSEAQSLKLGYVERWFLGSPLTVGFDFELTHKNLFVYRAGAKGNGLPHPYVSKEHWANSPGLAESFRLKYSRFESAIGAHTGYQWYPRYAVIRVNGGVDFRVVKNFYDKDNNQPFDLTVKEQLNWTSINSFWTSVSFDGRDFAYDPSSGWFLGQRCTFNGLVPCLEKEHSFRSDTKAEFYVTLLNYPVSAVWNLKFVLAFYTGVSVQTYYGRRKSENGKGNGVRSGALVIDGVLVGRGWSEDAKKNTGDLLLHHWIEFRWPLAHGIVSFDFFFDAAMVYNIESQSPNGSSSASSSSSSSSSSSSTTSS----EGLYKMSYGPGLRFTLPQFPLKLAFANTFTSPGGIPKTKKDWNFVLSFTVNNL";
        String alternative = "MLKKASAFLIASCCVMSLAWAQANDNWYEGKPISAISFEGLEYIARGQLDTIFSQYKGQKWTYELYLEILQKVYDLEYFSEVSPKAVPTDPEYQYVMLQFTVKERPSVKGIKMVGNSQIRSGDLLSKILLKKGDIYNEVKMKVDQESLRRHYLDQGYAAVKISCEAKTEAGGVVVQFTIQEGKQTVVSRIQFKGNKAFTESVLKKVLSTQEARFLTSGVFKENALEADKAAVHSYYAERGYIDARVEGVAKTVDKKTDASRNLVTLTYTVVEGEQYRYGGVTIVGNQIFSTEELQAKIRLKRGAIMNMVAFEQGFQALADAYFENGYTSNYLNKEEHRDTAEKTLSFKITVVERERSHVEHIIIKGTKNTKDEVILREMLLKPGDVFSKSKFTDSLRNLFNLRYFSSLVPDVRPGSEQDLVDIILNVEEQSTANVQFGVTFSGVGEAGTFPLSLFCQWEEKNFLGKGNEISVNATLGSEAQSLKLGYVERWFLGSPLTVGFDFELTHKNLFVYRAGSYGNGLPHPYTSREQWASSPGLAESFRLKYSRFESAIGAHTGYQWYPRYAVIRVNGGVDFRVVKNFYDKDNNQPFDLTVEEQLNWTSINSFWTSVSFDGRDFAYDPSSGWFLGQRCTFNGLVPFLEKEHSFRSDTKAEFYVTLLNYPVSAVWNLKFVLAFYTGVSVQTYYGRRKSENGKGNGVRSGALVIDGVLVGRGWSEDAKKNTGDLLLHHWIEFRWPLAHGIVSFDFFFDAAMVYNIESQSPNGSSSASSSSSSSSSSSSSSSSSSSSEGLYKMSYGPGLRFTLPQFPLKLAFANTFTSPGGIPKTKKNWNFVLSFTVNNL";
        ArrayList<Triple<Integer, String, String>> result = SequenceOperations.getCanonicalVariants(reference, alternative);
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
    void getCanonicalVariants_exception() {
        String reference = "GTATGGGGCT";
        String alternative = "GGGGGCT";
        assertThrowsExactly(IllegalArgumentException.class, () -> SequenceOperations.getCanonicalVariants(reference, alternative));
    }

    @Test
    void translateSequence_validSequence() throws MusialException {
        String sequence = "ATGCGT";
        String result = SequenceOperations.translateSequence(sequence, false);
        assertEquals("MR", result);
    }

    @Test
    void translateSequence_reverseSequence() throws MusialException {
        String sequence = "ATGCGT";
        String result = SequenceOperations.translateSequence(sequence, true);
        assertEquals("TH", result);
    }

    @Test
    void translateSequence_emptySequence() throws MusialException {
        String sequence = "";
        String result = SequenceOperations.translateSequence(sequence, false);
        assertEquals("", result);
    }

    @Test
    void translateSequence_invalidSequence() {
        String sequence = "ATGCGTX";
        assertThrowsExactly(MusialException.class, () -> SequenceOperations.translateSequence(sequence, false));
    }
}
