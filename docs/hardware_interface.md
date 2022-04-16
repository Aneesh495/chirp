# Hardware Interface and Sensor Telemetry Specification

## 1. Overview

CHIRP provides native ingest capabilities for Texas Instruments (TI) millimeter-wave radar front-ends paired with high-speed data capture cards:
- **RF Front-Ends**: TI AWR1243, AWR1843, AWR2243 (76 to 81 GHz automotive transceivers).
- **Data Capture Card**: TI DCA1000EVM (High-Speed Raw Data Capture FPGA Board).
- **Physical Link**: Gigabit Ethernet (1000BASE-T RJ45), UDP/IP streaming.

## 2. DCA1000EVM Packet Framing Protocol

The DCA1000EVM FPGA packetizes deserialized LVDS lane data and transmits raw baseband samples over Ethernet using User Datagram Protocol (UDP).

### UDP Packet Structure
Each Ethernet frame carries a fixed 1456-byte payload comprising a 10-byte hardware header followed by 1446 bytes of raw packed ADC samples:

| Byte Offset | Field Name | Data Type | Description |
| :--- | :--- | :--- | :--- |
| `0x00 .. 0x03` | `SequenceNumber` | `uint32` (little-endian) | Monotonically increasing packet sequence ID |
| `0x04 .. 0x09` | `ByteOffset` | `uint48` / `uint64` (little-endian)| Cumulative raw ADC bytes transmitted |
| `0x0A .. 0x5AF` | `RawSamplePayload`| `int16` array (1446 bytes) | Consecutive packed complex baseband I/Q samples |

### Network Configuration
- **Host IP**: `192.168.33.30`
- **DCA1000 IP**: `192.168.33.180`
- **Data Port**: `4098` (UDP stream from FPGA)
- **Config Port**: `4096` (Command and control JSON over UDP)

## 3. ADC Sample Demultiplexing and Channel Packing

TI mmWave sensors output 16-bit signed two's complement integers representing in-phase ($I$) and quadrature ($Q$) components.

### 4-Receiver Complex Format
For a 4-channel physical receive array, the FPGA transmits samples in continuous non-interleaved or interleaved blocks:

$$I_{\text{RX1}}[0], Q_{\text{RX1}}[0], I_{\text{RX2}}[0], Q_{\text{RX2}}[0], I_{\text{RX3}}[0], Q_{\text{RX3}}[0], I_{\text{RX4}}[0], Q_{\text{RX4}}[0], \dots$$

Each complex sample occupies 4 bytes (2 bytes $I$ + 2 bytes $Q$). A full frame requires:

$$\text{BytesPerFrame} = N_{\text{adc}} \times N_{\text{chirp}} \times N_{\text{RX}} \times 4\text{ bytes}$$

For a standard configuration ($N_{\text{adc}} = 256$, $N_{\text{chirp}} = 64$, $N_{\text{RX}} = 4$):

$$\text{BytesPerFrame} = 256 \times 64 \times 4 \times 4 = 262,144\text{ bytes} \approx 256\text{ KB}$$

At 20 frames per second, the sustained network throughput is $\approx 41.9\text{ Mbps}$, well within Gigabit Ethernet limits.

## 4. Packet Loss Detection and Stream Resynchronization

Because UDP is an unreliable transport protocol, high system CPU load or network buffering can cause dropped packets. The `Dca1000Reader` class implements strict packet continuity tracking:

```matlab
expectedSequence = lastSequence + 1;
if currentSequence ~= expectedSequence
    lostPackets = currentSequence - expectedSequence;
    warning('Dca1000Reader:PacketDrop', 'Detected %d dropped UDP packets.', lostPackets);
    % Zero-pad missing sample positions to preserve frame alignment
end
```

Frames with more than 1% lost packets are rejected to prevent corrupted Range-Doppler phase history.

## 5. Binary Radar Cube Container (`.rcube`) Specification

To store and exchange high-dimensional 3D and 4D radar cubes across operating systems without lossy floating-point conversion, CHIRP introduces the standardized binary `.rcube` container format:

### Header Format (64 Bytes)
| Byte Offset | Field | Type | Value / Description |
| :--- | :--- | :--- | :--- |
| `0x00 .. 0x03` | Magic Identifier | `char[4]` | `'RCUB'` (0x52, 0x43, 0x55, 0x42) |
| `0x04 .. 0x05` | Format Version | `uint16` | `0x0001` (Version 1.0) |
| `0x06 .. 0x07` | Compression Type | `uint16` | `0 = None`, `1 = Deflate`, `2 = Zstandard` |
| `0x08 .. 0x0B` | Number of ADC Samples ($N_R$) | `uint32` | Fast-time dimension |
| `0x0C .. 0x0F` | Number of Chirps ($N_V$) | `uint32` | Slow-time dimension |
| `0x10 .. 0x13` | Number of Virtual Channels ($N_C$)| `uint32` | Spatial co-array dimension |
| `0x14 .. 0x17` | Sample Data Type | `uint32` | `1 = Int16`, `2 = Float32`, `3 = Float64` |
| `0x18 .. 0x1F` | Center Frequency $f_c$ | `float64` | Carrier frequency in Hz |
| `0x20 .. 0x27` | Sweep Bandwidth $B$ | `float64` | Chirp bandwidth in Hz |
| `0x28 .. 0x2F` | Sampling Rate $f_s$ | `float64` | ADC sampling rate in Hz |
| `0x30 .. 0x33` | Payload CRC-32 | `uint32` | IEEE 802.3 32-bit cyclic redundancy check |
| `0x34 .. 0x3F` | Reserved Padding | `uint8[12]` | Zero-filled for 64-byte alignment |

The binary header is immediately followed by raw continuous complex double or single precision floats ordered as `[N_adc x N_chirps x N_channels]`. The container is read and written via `chirp.io.RadarCubeCodec`.
