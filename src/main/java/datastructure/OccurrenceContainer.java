package datastructure;

import java.util.Iterator;
import java.util.LinkedHashSet;

/**
 * Interface for sample occurrence storage.
 *
 * @author Simon Hackl
 * @version 2.3
 * @since 2.3
 */
public class OccurrenceContainer extends InfoContainer {

    /**
     * Set of all samples, represented by their internal name, allocated to this instance.
     */
    private final LinkedHashSet<String> occurrence = new LinkedHashSet<>(10, 10);

    /**
     * Constructor of {@link OccurrenceContainer}.
     */
    public OccurrenceContainer() {
        super();
    }

    /**
     * @param sampleName Adds specified name to this instance.
     */
    public void addOccurrence(String sampleName) {
        this.occurrence.add(sampleName);
    }

    /**
     * @param sampleName Sample name for which the occurrence of this instance is to be checked.
     * @return True, if {@code sampleName} is stored in this instance.
     */
    public boolean hasOccurrence(String sampleName) {
        return this.occurrence.contains(sampleName);
    }

    /**
     * @return Iterator over all stored values in this instance.
     */
    public Iterator<String> getOccurrenceIterator() {
        return this.occurrence.iterator();
    }

    /**
     * @return Set view of this instance.
     */
    public LinkedHashSet<String> getOccurrenceSet() {
        return this.occurrence;
    }

    /**
     * @return Size of this instance.
     */
    public int getOccurrenceCount() {
        return this.occurrence.size();
    }

}
