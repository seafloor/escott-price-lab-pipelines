skip_on_ci <- function() {
  if (Sys.getenv("CI") != "") {
    skip("Skipping test in CI environment")
  }
}