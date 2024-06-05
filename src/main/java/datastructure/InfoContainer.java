package datastructure;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Interface for meta-information storage.
 *
 * @author Simon Hackl
 */
@SuppressWarnings("unused")
public class InfoContainer {

    /**
     * Stores arbitrary meta information.
     */
    private final HashMap<String, String> info = new HashMap<>(5, 1);

    /**
     * Constructor of {@link InfoContainer}.
     */
    public InfoContainer() {
    }

    /**
     * Add meta information ({@code key}: {@code value}) to this instance.
     *
     * @param key   Key under which the value is stored.
     * @param value Value to be stored under Key.
     */
    public void addInfo(String key, String value) {
        this.info.put(key, value);
    }

    /**
     * Remove meta information stored at {@code key} from this instance.
     *
     * @param key Key whose value is to be deleted.
     */
    public void removeInfo(String key) {
        this.info.remove(key);
    }

    /**
     * Returns the value stored at {@code key} in this instance.
     *
     * @param key Key whose value is to be returned.
     * @return Value stored at {@code key} or {@code null} if no value is stored.
     */
    public String getInfo(String key) {
        return this.info.getOrDefault(key, null);
    }

    /**
     * @return All meta information stored for this instance.
     */
    public Set<Map.Entry<String, String>> getInfoSet() {
        return this.info.entrySet();
    }

    /**
     * @return All stored keys of this instance.
     */
    public Set<String> getInfoKeys() {
        return this.info.keySet();
    }

    /**
     * @param key Key for which it is to be checked whether it is stored.
     * @return True, if {@code key} is stored.
     */
    public boolean hasInfoKey(String key) {
        return this.info.containsKey(key);
    }

}
