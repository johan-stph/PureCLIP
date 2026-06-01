class Pureclip2 < Formula
  desc "Detect protein-RNA interaction footprints from single-nucleotide CLIP-seq data"
  homepage "https://github.com/johan-stph/PureCLIP"
  url "https://github.com/johan-stph/PureCLIP/archive/refs/tags/v2.0.0.tar.gz"
  sha256 "4583a4b6061264955c5d5194562979c6fa5f7301f224a7718eae66c69b2c0531"
  license "MIT"

  depends_on "cmake" => :build
  depends_on "boost"
  depends_on "gsl"
  depends_on "libomp"

  def install
    mkdir "build" do
      system "cmake", "../src",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DGSL_ROOT_DIR=#{Formula["gsl"].opt_prefix}",
        "-DCMAKE_PREFIX_PATH=#{Formula["boost"].opt_prefix}",
        *std_cmake_args
      system "make", "-j#{ENV.make_jobs}"
    end
    bin.install "build/pureclip"
    bin.install "build/winextract"
  end

  test do
    system "#{bin}/pureclip", "--version"
  end
end
