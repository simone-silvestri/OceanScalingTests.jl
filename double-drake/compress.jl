using OceanScalingTests
OceanScalingTests.compress_all_restarts((1080, 450, 100), 8, "./"; prefix = "DoubleDrake_checkpoint_")
OceanScalingTests.compress_surface_fields((1080, 450, 100), 8, "./"; prefix = "DoubleDrake_fields_")
