function floquet_drive(t::Float64, amplitude::Float64, omega::Float64, phase=0.0)
    return amplitude * cos(omega * t + phase)
end