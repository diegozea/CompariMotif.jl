using TestItemRunner

# Enforce oracle-independent tests regardless of local environment setup.
withenv("SLiMSuite_PATH" => nothing) do
    @run_package_tests
end
