package model;

import org.ehcache.spi.serialization.SerializerException;
import org.junit.jupiter.api.Test;

import java.nio.ByteBuffer;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class VariantCallSerializerTest {

    @Test
    void serializeAndDeserializeVariantCallShouldBeEqual() throws SerializerException {
        VariantCall original = new VariantCall(
                VariantCall.Flag.PASS,
                (short) 100,
                0.95f,
                List.of(new VariantCall.CallAlternative("A", "T", (short) 50))
        );

        VariantCall.VariantCallSerializer serializer = new VariantCall.VariantCallSerializer();
        ByteBuffer buffer = serializer.serialize(original);
        VariantCall deserialized = serializer.read(buffer);

        assertEquals(original, deserialized);
    }

    @Test
    void serializeWithEmptyAlternativesShouldWork() throws SerializerException {
        VariantCall original = new VariantCall(
                VariantCall.Flag.LOW_COVERAGE,
                (short) 0,
                0.0f,
                List.of()
        );

        VariantCall.VariantCallSerializer serializer = new VariantCall.VariantCallSerializer();
        ByteBuffer buffer = serializer.serialize(original);
        VariantCall deserialized = serializer.read(buffer);

        assertEquals(original, deserialized);
    }

    @Test
    void deserializeWithCorruptedDataShouldThrowException() {
        ByteBuffer corruptedBuffer = ByteBuffer.allocate(10);
        corruptedBuffer.putInt(5); // Invalid flag length
        corruptedBuffer.flip();

        VariantCall.VariantCallSerializer serializer = new VariantCall.VariantCallSerializer();

        assertThrows(SerializerException.class, () -> serializer.read(corruptedBuffer));
    }

    @Test
    void equalsShouldReturnTrueForMatchingVariantCallAndBinary() throws SerializerException {
        VariantCall original = new VariantCall(
                VariantCall.Flag.PASS,
                (short) 100,
                0.95f,
                List.of(new VariantCall.CallAlternative("A", "T", (short) 50))
        );

        VariantCall.VariantCallSerializer serializer = new VariantCall.VariantCallSerializer();
        ByteBuffer buffer = serializer.serialize(original);

        assertTrue(serializer.equals(original, buffer));
    }

    @Test
    void equalsShouldReturnFalseForNonMatchingVariantCallAndBinary() throws SerializerException {
        VariantCall original = new VariantCall(
                VariantCall.Flag.PASS,
                (short) 100,
                0.95f,
                List.of(new VariantCall.CallAlternative("A", "T", (short) 50))
        );

        VariantCall different = new VariantCall(
                VariantCall.Flag.LOW_FREQUENCY,
                (short) 50,
                0.5f,
                List.of(new VariantCall.CallAlternative("G", "C", (short) 25))
        );

        VariantCall.VariantCallSerializer serializer = new VariantCall.VariantCallSerializer();
        ByteBuffer buffer = serializer.serialize(different);

        assertFalse(serializer.equals(original, buffer));
    }
}
