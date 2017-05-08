# Coupled_gas_water_transport_with_geochemistry
A two-phase (water + gas) reactive transport model, written in Julia, that also simulates geochemical reactions between water-gas-mineral phases via a link with PHREEQC.

This code is an enhanced version of my gas_water_trans.jl script, https://github.com/NumericalEnvironmental/Coupled_gas_water_flow_and_transport_porous_media, with an interface with PHREEQC added to model equilibrium geochemical reactions between the water and gas phases in each cell, in the (possible) presence of a mineral phase assemblage. Julia version 0.5 or greater is required. In addition to some minor modifications to the gas_water_trans.jl module, this new code contains two additional modules:

* geochemical_module.jl - contains routines to write input for PHREEQC, run PHREEQC (via the IPhreeqcCOM library, https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/, which limits the model to Windows machines for the time being), and read PHREEQC's output

* header_phreeqc.jl - contains key word blocks to be repeated for all PHREEQC runs (all cells, all time steps); examples would include phase definitions

In addition to the input files required by the gas_water_trans.jl script, this code also requires the minerals.dat file, which defines the initial mineral assemblage for each node, to be present in the folder. Set the chemrxt boolean flag at the end of the knobs.txt file to enable the geochemical modeling (and be prepared for the model to run much slower ...).

This is just a test/demo prototype and has undergone only very limited testing. It may contain some bugs that could manifest in certain complex or large reactive transport problems. It will also likely be very slow for larger problems; I definitely need to do some more work with time stepping. A little more discussion is provided in my blog, (link coming shortly).

Again, please note that the IPhreeqcCOM library MUST be present on the Windows machine where this script is run. In addition, the PyCall package for Julia must also be installed.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
