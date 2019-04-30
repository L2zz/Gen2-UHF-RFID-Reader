"""
Decode EPC data with several algorithms and get result
"""
import sys
import csv
import time

"""
TODO
: Get data start from last half bit of preamble
: Testing
"""

SLIDING_WINDOW = 3
NUM_OF_ROUNDS = 1000
NUM_OF_SAMPLES_PER_ONE_BIT = 49
NUM_OF_EPC_BITS = 128

class DecodingAlgorithm:
    """
    Define decoding algorithms
    """
    NUM_OF_ALGORITHM = 3
    INSPECT_ALL = 0
    PRESIZE = 1
    THRESHOLD = 2
    BIT_MASK = [[[-1, 1, -1, 1], [-1, 1, 1, -1]],
                [[1, -1, 1, -1], [1, -1, -1, 1]]]
    DECODING_THRESHOLD = 0.8

    @classmethod
    def get_first_masklevel(cls, data, single_unit):
        """
        Determine mask level of first bit
        """
        mask_level = 1
        max_corr = 0.0
        for i in range(2):
            index = single_unit
            _, corr = cls.decoding_single_bit(data, index, i, single_unit)
            if corr > max_corr:
                max_corr = corr
                mask_level = i
        if mask_level == 0:
            mask_level = -1

        return mask_level

    @classmethod
    def get_presize(cls, data, index):
        """
        Determine size of full signal
        """
        avg_amp = 0.0
        expected_end = NUM_OF_SAMPLES_PER_ONE_BIT*NUM_OF_EPC_BITS
        for i in range(int(-3.5*NUM_OF_SAMPLES_PER_ONE_BIT),
                       int(0.5*NUM_OF_SAMPLES_PER_ONE_BIT)):
            avg_amp += data[expected_end+i]
        avg_amp /= 4*NUM_OF_SAMPLES_PER_ONE_BIT

        end = expected_end
        for i in range(int(-3.5*NUM_OF_SAMPLES_PER_ONE_BIT),
                       int(0.5*NUM_OF_SAMPLES_PER_ONE_BIT)):
            prev_amp = data[expected_end+i-1]
            cur_amp = data[expected_end+i]
            if (cur_amp-avg_amp)*(prev_amp-avg_amp) < 0:
                end = expected_end + i
        presize = (end-index) / (NUM_OF_EPC_BITS*1.0)

        return presize

    @classmethod
    def decoding_single_bit(cls, data, index, prev_level, single_unit):
        """
        Determine single bit using correlation comparison
        """
        avg_amp = 0.0
        for i in range(int(-0.5*single_unit), int(1.5*single_unit)):
            avg_amp += abs(data[index+i])
        avg_amp /= 2.0*single_unit
        avg_abs_amp = 0.0
        for i in range(int(-0.5*single_unit), int(1.5*single_unit)):
            avg_abs_amp += abs(data[index+i] - avg_amp)
        avg_abs_amp /= 2.0*single_unit

        single_bit_decoded = 0
        max_corr = 0.0
        if prev_level == -1:
            prev_level = 0
        for i in range(2):
            corr = 0.0
            for j in range(int(-0.5*single_unit), int(1.5*single_unit)):
                if j < 0:
                    half_bit_idx = 0
                elif j < int(0.5*single_unit):
                    half_bit_idx = 1
                elif j < int(single_unit):
                    half_bit_idx = 2
                else:
                    half_bit_idx = 3
                scaled_amp = (data[index+j] - avg_amp) / avg_abs_amp
                corr += cls.BIT_MASK[prev_level][i][half_bit_idx] * scaled_amp
            if corr > max_corr:
                max_corr = corr
                single_bit_decoded = i
        max_corr /= 2.0*single_unit

        return single_bit_decoded, max_corr

    @classmethod
    def decoding_inspect_all(cls, data):
        """
        Decoding scheme by inspecting all the candidates
        """
        one_round_decoded = []
        shift = 0
        prev_level = cls.get_first_masklevel(data, NUM_OF_SAMPLES_PER_ONE_BIT)
        for i in range(NUM_OF_EPC_BITS):
            max_corr = 0.0
            index = i*NUM_OF_SAMPLES_PER_ONE_BIT + shift - SLIDING_WINDOW
            for j in range(SLIDING_WINDOW*2+1):
                index += j
                bit_candidate, corr = cls.decoding_single_bit(data, index, prev_level,
                                                              NUM_OF_SAMPLES_PER_ONE_BIT)
                if corr > max_corr:
                    max_corr = corr
                    decoded_bit = bit_candidate
                    cur_shift = j - SLIDING_WINDOW
            one_round_decoded.append(decoded_bit)
            if decoded_bit == 1:
                prev_level *= -1
            shift += cur_shift

        return one_round_decoded

    @classmethod
    def decoding_presize(cls, data):
        """
        Decoding scheme using empirical size of signal
        """
        one_round_decoded = []
        presize_samples = cls.get_presize(data, 0)
        prev_level = cls.get_first_masklevel(data, presize_samples)
        for i in range(NUM_OF_EPC_BITS):
            index = i*presize_samples
            decoded_bit, _ = cls.decoding_single_bit(data, index, prev_level,
                                                     presize_samples)
            one_round_decoded.append(decoded_bit)
            if decoded_bit == 1:
                prev_level *= -1

        return one_round_decoded

    @classmethod
    def decoding_threshold(cls, data):
        """
        Decoding scheme using threshold
        """
        one_round_decoded = []
        shift = 0
        prev_level = cls.get_first_masklevel(data, NUM_OF_SAMPLES_PER_ONE_BIT)
        for i in range(NUM_OF_EPC_BITS):
            max_corr = 0.0
            cur_shift = 0
            index = i*NUM_OF_SAMPLES_PER_ONE_BIT + shift - SLIDING_WINDOW
            decoded_bit, max_corr = cls.decoding_single_bit(data, index, prev_level,
                                                            NUM_OF_SAMPLES_PER_ONE_BIT)
            if max_corr < cls.DECODING_THRESHOLD:
                for j in range(SLIDING_WINDOW*2+1):
                    index += j
                    bit_candidate, corr = cls.decoding_single_bit(data, index, prev_level,
                                                                  NUM_OF_SAMPLES_PER_ONE_BIT)
                    if corr > max_corr:
                        max_corr = corr
                        decoded_bit = bit_candidate
                        cur_shift = j - SLIDING_WINDOW
            one_round_decoded.append(decoded_bit)
            if decoded_bit == 1:
                prev_level *= -1
            shift += cur_shift

        return one_round_decoded

def print_result(source_file, decode_result, time_avg_result):
    """
    Print results
    """
    result_file = open("result", 'a')
    csv_wr = csv.writer(result_file, delimiter=',')

    result_arr = []
    result_arr.append(source_file)

    total = 0
    success = 0
    for result in decode_result:
        if result:
            success += 1
        total += 1
    accuracy = (success*1.0) / total
    result_arr.append(accuracy)

    time_avg = 0.0
    for time_result in time_avg_result:
        time_avg += time_result
    time_avg /= (1.0*len(time_avg_result))

    result_arr.append(time_avg)
    csv_wr.writerow(result_arr)

    print(result_arr + ": Accuracy = %.4f Time avg = %.4f" % (accuracy, time_avg))
    result_file.close()

def check_crc(decoded_data):
    """
    CRC checking for EPC data
    """
    char_bits = []
    for i in range(NUM_OF_EPC_BITS):
        if decoded_data[i] == 0:
            char_bits.append('0')
        else:
            char_bits.append('1')

    char_data = [0 for _ in range(NUM_OF_EPC_BITS/8)]
    for i in range(NUM_OF_EPC_BITS/8):
        mask = 0x80
        char_data[i] = 0
        for j in range(8):
            if char_bits[8*i+j] == '1':
                char_data[i] = char_data[i] | mask
            mask = mask >> 1
    rcvd_crc = (char_data[NUM_OF_EPC_BITS/8 - 2] << 8) + char_data[NUM_OF_EPC_BITS/8 - 1]

    crc_16 = 0xFFFF
    for i in range(NUM_OF_EPC_BITS/8 - 2):
        crc_16 = crc_16 ^ (char_data[i] << 8)
        for j in range(8):
            if crc_16 & 0x8000:
                crc_16 = crc_16 << 1
                crc_16 = crc_16 ^ 0x1021
            else:
                crc_16 = crc_16 << 1
    crc_16 = ~crc_16

    if rcvd_crc != crc_16:
        result = False
    else:
        result = True

    return result

def decode_data(data, decoding_id):
    """
    Decode data
    param:
        data: Data array to decode
        decoding_id: Unique id of the decoding algorithm to use
    return:
        decoded_data: Decoding result array
        decoding_time_avg: Average execution time of the target decoding algorithm
    """
    decoded_data = []
    decoding_time_avg = 0.0
    if decoding_id == DecodingAlgorithm.INSPECT_ALL:
        for i in range(NUM_OF_ROUNDS):
            start_time = time.time()
            one_round_decoded = DecodingAlgorithm.decoding_inspect_all(data[i])
            decoded_data.append(one_round_decoded)
            decoding_time_avg += time.time() - start_time
        decoding_time_avg /= NUM_OF_ROUNDS*1.0
    elif decoding_id == DecodingAlgorithm.PRESIZE:
        for i in range(NUM_OF_ROUNDS):
            start_time = time.time()
            one_round_decoded = DecodingAlgorithm.decoding_presize(data[i])
            decoded_data.append(one_round_decoded)
            decoding_time_avg += time.time() - start_time
        decoding_time_avg /= NUM_OF_ROUNDS*1.0
    elif decoding_id == DecodingAlgorithm.THRESHOLD:
        for i in range(NUM_OF_ROUNDS):
            start_time = time.time()
            one_round_decoded = DecodingAlgorithm.decoding_threshold(data[i])
            decoded_data.append(one_round_decoded)
            decoding_time_avg += time.time() - start_time
        decoding_time_avg /= NUM_OF_ROUNDS*1.0
    else:
        print("Unknown algorithm")
        sys.exit()

    return decoded_data, decoding_time_avg

def get_data(source_file):
    """
    Get target data from source_file
    param:
        source_file: File which includes data written with text form
    return:
        data: Data array
    """
    src = open(source_file, 'r')
    data = []
    for _ in range(NUM_OF_ROUNDS):
        one_round = src.readline()
        if not one_round:
            print("Not enough data")
            sys.exit()
        one_round_data = one_round.split(' ')
        one_round_data.pop()
        data.append(one_round_data)
    src.close()

    return data

def get_result(source_files):
    """
    Get results of decoding with several algorithms
    param:
        source_files: File array to read
    """
    for source_file in source_files:
        time_avg_arr = []
        decoding_result = []
        data = get_data(source_file)
        for i in range(DecodingAlgorithm.NUM_OF_ALGORITHM):
            decoded_data, time_avg = decode_data(data, i)
            decoding_result.append(check_crc(decoded_data))
            time_avg_arr.append(time_avg)
        print_result(source_file, decoding_result, time_avg_arr)

if __name__ == "__main__":

    SOURCE_FILES = ["100_0_0"]
    get_result(SOURCE_FILES)
