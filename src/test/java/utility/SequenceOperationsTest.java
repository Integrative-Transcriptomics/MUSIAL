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
        String sequenceA = "ACGT";
        String sequenceB = "ACGT";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACGT", "ACGT"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleSubstitution() {
        String sequenceA = "ACGT";
        String sequenceB = "AGGT";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACGT", "AGGT"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleInsertion() {
        String sequenceA = "ACGT";
        String sequenceB = "ACGTT";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACGT-", "ACGTT"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_singleDeletion() {
        String sequenceA = "ACGT";
        String sequenceB = "ACT";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACGT", "AC-T"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_emptySequenceA() {
        String sequenceA = "";
        String sequenceB = "ACGT";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("----", "ACGT"), result);
    }

    @Test
    void globalNucleotideSequenceAlignment_emptySequenceB() {
        String sequenceA = "ACGT";
        String sequenceB = "";
        Tuple<String, String> result = SequenceOperations.globalNucleotideSequenceAlignment(sequenceA, sequenceB, 2, 1,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACGT", "----"), result);
    }

    @Test
    void globalProteinSequenceAlignment_identicalSequences() {
        String sequenceA = "ACDEFGHIKLMNPQRSTVWY";
        String sequenceB = "ACDEFGHIKLMNPQRSTVWY";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, 0);
        assertEquals(new Tuple<>("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"), result);
    }

    @Test
    void globalProteinSequenceAlignment_singleSubstitution() {
        String sequenceA = "ACDEFGHIKLMNPQRSTVWY";
        String sequenceB = "ACDEFGHIKLMNPQRSTVWZ";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, 0);
        assertEquals(new Tuple<>("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWZ"), result);
    }

    @Test
    void globalProteinSequenceAlignment_insertion() {
        String sequenceA = "ACDEFGHIKLMNPQRSTVWY";
        String sequenceB = "ACDEFGSSSHIKLMNPQRSTVWY";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, 3);
        assertEquals(new Tuple<>("ACDEFG---HIKLMNPQRSTVWY", "ACDEFGSSSHIKLMNPQRSTVWY"), result);
    }

    @Test
    void globalProteinSequenceAlignment_deletion() {
        String sequenceA = "ACDEFGHIKLMNPQRSTVW";
        String sequenceB = "ACDEFGHIKLQRSTVW";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, 3);
        assertEquals(new Tuple<>("ACDEFGHIKLMNPQRSTVW", "ACDEFGHIKL---QRSTVW"), result);
    }

    @Test
    void globalProteinSequenceAlignment_indel() {
        String sequenceA = "ACDEFGSSSHIKLMNPQRSW";
        String sequenceB = "ACDEFGHIKLMNPQRSTVW";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, 5);
        assertEquals(new Tuple<>("ACDEFGSSSHIKLMNPQRS--W", "ACDEFG---HIKLMNPQRSTVW"), result);
    }

    @Test
    void globalProteinSequenceAlignment_emptySequenceA() {
        String sequenceA = "";
        String sequenceB = "ACDEFGHIKLMNPQRSTVWY";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("--------------------", "ACDEFGHIKLMNPQRSTVWY"), result);
    }

    @Test
    void globalProteinSequenceAlignment_emptySequenceB() {
        String sequenceA = "ACDEFGHIKLMNPQRSTVWY";
        String sequenceB = "";
        Tuple<String, String> result = SequenceOperations.globalProteinSequenceAlignment(sequenceA, sequenceB, 8, 7,
                SequenceOperations.MarginalGaps.FORBID, SequenceOperations.MarginalGaps.PENALIZE, null);
        assertEquals(new Tuple<>("ACDEFGHIKLMNPQRSTVWY", "--------------------"), result);
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
