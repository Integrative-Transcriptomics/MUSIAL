package utilities;
import java.util.Objects;

public class Pair<F extends Comparable<F>, S extends Comparable<S>>
  implements Comparable<Pair<F, S>> {

    public final F first;
    public final S second;

    /**
     * Constructor for a Pair.
     *
     * @param first  the first object in the Pair
     * @param second the second object in the pair
     */
    public Pair(F first, S second) {
        this.first = first;
        this.second = second;
    }

    /**
     * Checks the two objects for equality by delegating to their respective
     * {@link Object#equals(Object)} methods.
     *
     * @param o the {@link Pair} to which this one is to be checked for equality
     * @return true if the underlying objects of the Pair are both considered
     * equal
     */
    @Override
    public boolean equals(Object o) {
        if (!(o instanceof Pair)) {
            return false;
        }
        Pair<?, ?> p = (Pair<?, ?>) o;
        return Objects.equals(p.first, first) && Objects.equals(p.second, second);
    }

    /**
     * Compute a hash code using the hash codes of the underlying objects
     *
     * @return a hashcode of the Pair
     */
    @Override
    public int hashCode() {
        return (first == null ? 0 : first.hashCode()) ^ (second == null ? 0 : second.hashCode());
    }

    /**
     * Convenience method for creating an appropriately typed pair.
     *
     * @param a the first object in the Pair
     * @param b the second object in the pair
     * @return a Pair that is templatized with the types of a and b
     */
    public static <A extends Comparable<A>, B extends Comparable<B>> Pair<A, B> create(A a, B b) {
        return new Pair<>(a, b);
    }

    @Override
    public int compareTo(Pair<F, S> that) {
        int cmp = this.first.compareTo(that.first);
        if (cmp == 0)
            cmp = this.second.compareTo(that.second);
        return cmp;
    }



	public F getFirst(){
		return this.first;
	}

	public S getSecond(){
		return this.second;
	}

	@Override
	public String toString(){
		return "("+this.first+","+this.second+")";
	}
}