
import os
import wave
import glob


for file_name in glob.glob('../*.wav'):
    with wave.open(file_name, "rb") as wave_file:
        frame_rate = wave_file.getframerate()
        print(file_name)
        print("Sample Rate %s" % frame_rate)


