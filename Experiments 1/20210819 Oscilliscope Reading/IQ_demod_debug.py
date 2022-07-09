if __name__ == '__main__':
    from pathlib import Path
    import numpy as np
    import matplotlib.pyplot as plt
    from IQ_demod import *

    def test_case_0():
        # checking if IQ demod for generalised case works
        SIGNAL_F = 25e3 #MHz
        SAMPLING_F = 100e3 #MHz

        sampling_interval = 1/SAMPLING_F
        print(f"{sampling_interval =}")
        signal = np.fromiter([np.sin(i * sampling_interval) for i in range(150)], dtype=np.float64)
        
        plt.plot(signal)
        plt.show(block=True)
        plt.pause(0.8)

        N = 4 # number of initialisation points wanted
        delta_phi = phase_advance(SIGNAL_F, SAMPLING_F)
        print(f"{delta_phi = }")
        
        phases = signal_to_phase(signal, N, delta_phi, True)

        # phase reconstruction 1
        phases_reconstructed = [phases[0]]
        mod_twopi = 0
        for i in range(1, len(phases)):
            if phases[i] - phases[i-1] > 4.5:
                # phases_reconstructed.append(phases_reconstructed[i-1] + phases[i] - 2*pi)
                mod_twopi -= 1
            elif phases[i] - phases[i-1] < -4.5:
                mod_twopi += 1
            phases_reconstructed.append(phases[i] + 2 * pi * mod_twopi)
        plt.plot(phases_reconstructed)
        plt.draw()
        plt.pause(0.8)

    def test_case_1():
        # attempting to do from scratch
        # TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210812 Oscilliscope Matters\\Test Files\\C2Test16Aug202100001.txt')
        # TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Reference\\100ks\\C2-162mhz-100ks-00000.txt')
        TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Reference\\C2-ball-drop-00002.txt')
        readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        meta, trace = fr.parse_oscilliscope_txt(readTxt)
        print('meta = ')
        pprint(meta)
        print()

        signal = fr.signal_from_trace(trace)
        signal = signal[:100000]
        # signal = np.array(signal)

        # plt.plot(signal)
        # plt.show(block=False)

        # frequencies
        SIGNAL_F = 160.00025e6 #MHz
        SAMPLING_F = 100e3 #MHz

        # trace_intervals = int(meta['Sample Interval'][0] * SAMPLING_F)
        # print(f"{trace_intervals = }")
        print(f"sampling interval = {1/SAMPLING_F}")

        # # Initialisation parameters
        N = 4 # number of initialisation points wanted
        delta_phi = phase_advance(SIGNAL_F, SAMPLING_F)
        print(f"{delta_phi = }")
        
        phases = signal_to_phase(signal, N, delta_phi, True)
        # tolerance for when to wrap over phases
        tol = 4.8

        # phase reconstruction 1
        phases_reconstructed = [phases[0]]
        mod_twopi = 0
        for i in range(1, len(phases)):
            if phases[i] - phases[i-1] > tol:
                # phases_reconstructed.append(phases_reconstructed[i-1] + phases[i] - 2*pi)
                mod_twopi -= 1
            elif phases[i] - phases[i-1] < -tol:
                mod_twopi += 1
            phases_reconstructed.append(phases[i] + 2 * pi * mod_twopi)
        plt.plot(phases_reconstructed)
        # plt.show(block=True)
        # plt.pause(0.8)

        phases_r = phase_reconstruction_naive(phases, tol)
        print(f"Debug: {phases[:15] = }")
        # print(f"Debug: {phases_reconstructed[:15] = }")

        plt.plot(phases_r)
        plt.show(block=True)
        plt.pause(0.8)

    def test_case_2a():
        # checking that the matlab version of IQ demod works as intended.
        TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Reference\\C2-ball-drop-00002.txt')
        readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        meta, trace = fr.parse_oscilliscope_txt(readTxt)
        print('meta = ')
        pprint(meta)
        print()

        signal = fr.signal_from_trace(trace)
        phases = IQ_signal_to_phase(signal)
        phases = phase_reconstruction(phases, 4.8)
        plt.plot(phases)
        plt.show(block = True)
        # the real version works faster than the naive version
    
    def test_case_2b():
        # checking that the matlab version of IQ demod works as intended.
        TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Reference\\C2-ball-drop-00002.txt')
        readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        meta, trace = fr.parse_oscilliscope_txt(readTxt)
        print('meta = ')
        pprint(meta)
        print()

        signal = fr.signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, 4, pi/2, True)
        phases = phase_reconstruction(phases, 5)
        plt.plot(phases)
        plt.show(block = True)
        
    def test_case_3():
        N, M = freq_ratio(160, 15)
        print(f"{N = }")
        print(f"{M = }")

    def test_case_4():
        # determine if can M/N appropriately directly from function 
        # freq_ratio
        # properly shows table knocking over time

        TXT_FILE_PATH = Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Reference\\C2-ball-drop-00002.txt')
        readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        meta, trace = fr.parse_oscilliscope_txt(readTxt)

        N, M= freq_ratio(signal=160.025e6, sample=100e3)

        signal = fr.signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, N, 2*pi/N, True)
        phases = phase_reconstruction(phases, 5)
        t_axis = np.arange(start= 0, 
            stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(t_axis, phases)

        ax.set_ylabel(r'$\phi$/rad', useTex = True)
        ax.set_xlabel(r'$t$/s', useTex = True)

        plt.show(block = True)

    def test_case_5():
        # checking if the generalized IQ works for N != 4

        signal_freq = 80.02e6 #Hz
        sample_freq = 100e3 #Hz     
        N, M = freq_ratio(signal=signal_freq, sample=sample_freq)
        print(f"[Debug] {N = }")

        TXT_FILE_PATH = Path("D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Nis5data\\C280.020Mhz100kS-nothing00001.txt")
        readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        meta, trace = fr.parse_oscilliscope_txt(readTxt)
        signal = fr.signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, N, 2*pi/N, True)
        phases = phase_reconstruction(phases, 5)
        plt.plot(phases)
        plt.show(block = True)

    def test_case_6():
        # determine as to what the issue might be with reading for when 
        initial_phase = 0.248*pi
        linear_phase = 0

        signal_freq = 80.0250e6 #Hz
        sample_freq = 100e3 #Hz   
        N, M = freq_ratio(signal=signal_freq, sample=sample_freq)
        print(f"[Debug] {N = }")

        time_interval = 1/sample_freq
        iters = int(1 // time_interval)+1

        signal = np.sin(signal_freq * np.arange(1, iters) * time_interval + initial_phase)
        phases = signal_to_phase(signal, N, 2*pi/N, False)
        phases = phase_reconstruction(phases, 7)
        expected_phases = np.mod(np.arange(1, iters) * time_interval / (1/signal_freq) * 2*pi + initial_phase, 2*pi)
        expected_phases = phase_reconstruction(expected_phases, 7)
        plt.plot(phases)
        plt.plot(expected_phases)
        plt.show(block = True)
        # no real conclusive results

    def test_case_7():
        initial_phase = 0.138*pi
        quad_phase = -40

        signal_freq = 400 #Hz
        sample_freq = 300 #Hz   

        N, M = freq_ratio(signal=signal_freq, sample=sample_freq)
        print(f"[Debug] {N = }")

        k = 2.5
        signal_linspace = np.linspace(1,1+1/k, int(80*signal_freq/k))
        sample_linspace = np.arange(1,1+1/k,1/sample_freq)

        sample_signal = np.sin(signal_freq * sample_linspace + initial_phase + quad_phase * sample_linspace**2)
        signal_signal = np.sin(signal_freq * signal_linspace + initial_phase + quad_phase * signal_linspace**2)
        sample_phases = np.mod(signal_freq * sample_linspace + initial_phase + quad_phase * sample_linspace**2, 2*pi)
        signal_phases = np.mod(signal_freq * signal_linspace + initial_phase + quad_phase * signal_linspace**2, 2*pi)
        IQ_sample_phase = signal_to_phase(sample_signal, N, 2*pi/N, False)
        print(f"[Debug] {sample_signal[0:N] = }")
        IQ_sample_phase = IQ_sample_phase - (IQ_sample_phase[0] - signal_phases[0])
        IQ_sample_phase = np.asarray(list(map(lambda x: x if x >= 0 else 2*pi+x, IQ_sample_phase)))

        # waveform
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        ax1.scatter(sample_linspace, sample_signal, s = 20, c='g', marker = 'o')
        ax1.scatter(signal_linspace, signal_signal, s = 1, c = 'b')
        
        # phases
        ax2 = plt.subplot(212, sharex=ax1)
        ax2.scatter(signal_linspace, signal_phases/pi, s = 1, c = 'b')
        ax2.scatter(sample_linspace, sample_phases/pi, s = 15, c = 'g')
        ax2.scatter(sample_linspace[:len(sample_linspace)-N],IQ_sample_phase/pi, s = 35, c = 'r', marker = '1')
        

        fig.set_size_inches(11.75-1.5, 8.25-2 -0.6) # leave some space for table and section title
        
        ax1.set_ylabel(r'$y(t)$', useTex = True)
        ax2.set_ylabel(r'$\phi/\pi$ [rad]', useTex = True)
        ax2.set_xlabel(r'$t$ [s]', useTex = True)
        fig.suptitle("$f_{\\rm{signal}}$ = " + str(signal_freq) + 
            ", $f_{\\rm{sample}}$ = " + str(sample_freq) + ", $M/N$ = " 
            + str(M) + '/' + str(N), useTex = True)
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.95)
        plt.show(block = True)
        
    test_case_2b()
