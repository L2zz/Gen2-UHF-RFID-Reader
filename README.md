# Gen2 UHF RFID Reader with WISP
This program is revised version of the Gen2 UHF RFID Reader (https://github.com/nkargas/Gen2-UHF-RFID-Reader) developed by Nikos Kargas (https://github.com/nkargas). It can communicate with the programmable RFID Tag with FM0 encoding and 40kHz data rate.

## Hardware
 * 1x USRP N210
 * 1x SFX daughterboard
 * 2x circular polarized antennas

### Wireless Identification and Sensing Platform (WISP)
We use WISP(MSP430) as an RFID tag. Normal data rate of WISP is 460kHz, but we modified to 40kHz in order to communicate with this program. Visit below and install the 40kHz version on your WISP tag.
 * https://github.com/whitecloudy/wisp5_436kHz-40kHz/tree/40kHz

Get more information at below WISP wiki website.
 * http://wisp5.wispsensor.net

## Installation
1. Install the UHD driver and GNU Radio program. Use below instruction to check the success of the hardware installation.
<pre><code>$ uhd_find_devices`</code></pre>

2. Visit below and check the latest release of this program. Copy the link of the source code.
 * https://github.com/SpiritFlag/Gen2-UHF-RFID-Reader-with-WISP/releases

3. Download the source code.
<pre><code>$ wget https://github.com/SpiritFlag/Gen2-UHF-RFID-Reader-with-WISP/archive/(version).zip
$ wget https://github.com/SpiritFlag/Gen2-UHF-RFID-Reader-with-WISP/archive/(version).tar.gz</code></pre>

4. Unzip the ".zip" or ".tar.gz" file.
<pre><code>$ unzip *
$ tar -xvzf *</code></pre>

5. Execute the script file "init.sh" for initial build.
<pre><code>$ ./init.sh</code></pre>  
or type the below instructions manually.  
<pre><code>$ cd gr-rfid
$ mkdir build
$ cd build
$ cmake ../
$ make
$ make test
$ sudo make install
$ sudo ldconfig
$ cd ../misc
$ mkdir data
$ cd data
$ touch source matched_filter gate decoder
$ cd ../../apps</code></pre>

## Configuration
There are several variables you can modify in order to fit your experimental environment.

### gr-rfid/apps/reader.py
 * line 13: DEBUG  
True: execute the program from file source. (debugging mode)  
False : execute the program from USRP reader. (default)

 * line 53~59:  
dac_rate: DAC rate (default: 1MS/s)  
adc_rate: ADC rate (default: 2MS/s)  
decim: downsampling factor (default: 1)  
ampl: output signal amplitude (default: 0.55)  
freq: modulation frequency (default: 910MHz)  
rx_gain: RX gain  
tx_gain: TX gain

 * line 61~62:  
Change addr value with the address of your USRP reader. (default: 192.168.1.2)

### gr-rfid/include/global_vars.h
 * line 29: DEBUG_MESSAGE  
Annotate this line, if you don't want to make the debug files. (dafault: annotated)

 * line 75: FIXED_Q  
The number of slot is fixed by 2^(FIXED_Q). (default FIXED_Q: 0 / default slot number: 1)

 * line 79: MAX_NUM_QUERIES  
The program stops after sending this amount of queries. (dafault: 100000)

 * line 159: result_file_path  
Set the name of the result file. You have to change reader.sh file also. (default: result)

 * line 161: debug_message  
Set the name of the folder which saves the debug files. You have to change reader.sh file also. (dafault: debug_data/)

## Execution
Execute the "gr-rfid/apps/reader.py" python file. You must delete the "debug_data" folder before the every execution, because the program does not automatically remove the debug files from the previous execution. For convenience, there is a script file which automatically delete the unnecessary files. Use "reader.sh" rather than directly executing "reader.py".
<pre><code>$ ./reader.sh</code></pre>

If you made any change in ".cc" or ".h" file, you must rebuild the program. Execute the script file "build.sh". You don't need to rebuild the program when you only modify the "reader.py" file.
<pre><code>$ ./build.sh</code></pre>

or type the below instructions manually.
<pre><code>$ cd ../build && make && make test && sudo make install && sudo ldconfig && cd ../apps</code></pre>

When debugging mode, rename the "gr-rfid/misc/data/source" file to "file_source". Don't forget to make new "source" file after renaming.
<pre><code>$ touch source</code></pre>

If you want to reenact the execution, backup the "source" file in somewhere. You can easily reenact the execution by renaming file name to "file_source".

## Output
As the result of the execution, below files are created.

### Text File
 * gr-rfid/apps/debug_message  
Logs the flow of the program. It includes decoded RN16 bits and EPC bits.
 * gr-rfid/apps/result  
Logs the result of the program. It includes the detected tag IDs and the number of reads.
 * gr-rfid/apps/debug_data/(inventory_round)_(slot_number)  
Logs the squared normalized value of samples from each inventory round and slot number. It includes detailed flow of decoding process.
 * gr-rfid/apps/debug_data/(inventory_round)_(slot_number)_iq  
Logs the real and imaginary value of samples from each inventory round and slot number.

### Plot File
 * gr-rfid/misc/data/source  
Logs the all received samples. You should backup this file in order to reenact the execution.
 * gr-rfid/misc/data/matched_filter  
Logs the processed received samples. These samples are theinput of the reader block.
 * gr-rfid/misc/data/gate  
Logs the output of the gate block. These samples are the input of the tag_decoder block.
 * gr-rfid/misc/data/file_sink  
Logs the transmitted samples.

These files can be plotted by graphic interface. Use below instruction. The blue line figures the real value, and the red line figures the imaginary value.
<pre><code>$ gr_plot_iq -B (sample_window) (file_name)</code></pre>
For example, below instruction opens the matched_filter file with sample window 100000.
<pre><code>$ gr_plot_iq -B 100000 matched_filter</code></pre>

## Tested on:
Ubuntu 16.04 64-bit  
GNU Radio 3.7.10.1

## If you use this software please cite:
N. Kargas, F. Mavromatis and A. Bletsas, "Fully-Coherent Reader with Commodity SDR for Gen2 FM0 and Computational RFID", IEEE Wireless Communications Letters (WCL), Vol. 4, No. 6, pp. 617-620, Dec. 2015.

## Contacts
### Revised Version
 * Developed by  
Jeong SinGi (e-mail: jsgnwk@csi.skku.edu)  
Computer Science and Intelligence Lab in Sungkyunkwan University, Suwon, Republic of Korea(ROK, South Korea).

 * Co-worked with  
Shin JaeMin (e-mail: alex9395@csi.skku.edu)  
Computer Science and Intelligence Lab in Sungkyunkwan University, Suwon, Republic of Korea(ROK, South Korea).

 * Supervised by  
Prof. Kim YuSung (e-mail: yskim525@csi.skku.edu)  
Computer Science and Intelligence Lab in Sungkyunkwan University, Suwon, Republic of Korea(ROK, South Korea).

### Original Version
Nikos Kargas (e-mail1: cpznick@gmail.com e-mail2: karga005@umn.edu)
