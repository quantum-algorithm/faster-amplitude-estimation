from faster_amplitude_estimation import FasterAmplitudeEstimation, GroverSimulator

if __name__ == '__main__':
  estimation = FasterAmplitudeEstimation()
  amplitude = 0.1
  l = 8
  estimate, oracle_call, j_0 = estimation.estimate_amplitude(GroverSimulator(amplitude), l, 0.01)
  error = abs(estimate - amplitude)
  print(error)
