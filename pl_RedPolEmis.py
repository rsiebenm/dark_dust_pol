import numpy as np
import matplotlib.pyplot as plt

# Helper function to mimic IDL's rdfloat/readcol
def read_data(filename, skiprows=3):
    return np.loadtxt(filename, skiprows=skiprows, unpack=True)

# 1. Reddening Plot
w, sac, ssac, sasi, sssi, savgr, ssvgr, savsi, ssvsi, sapah, sad, ssd, st = read_data('./Output/Kappa.out')
w *= 1e4  # Convert to microns

kvis = np.where(w <= 0.551)[0][0]
kblue = np.where(w <= 0.45)[0][0]

av_mod = st[kvis] * 2.5 / np.log(10.0)
rv_mod = st[kvis] / (st[kblue] - st[kvis])
print(f'HD027778: Av_mod={av_mod:.3f}, Rv_mod={rv_mod:.3f}')

fig, axs = plt.subplots(2, 2, figsize=(14, 10))
plt.subplots_adjust(hspace=0.3, wspace=0.3)

# Plot 1: Reddening
ax1 = axs[0, 0]
inv_w = 1.0 / w
extinction = (st / st[kvis] - 1.0) * av_mod

ax1.fill_between(inv_w, -2, (sad + ssd) / st[kvis] * av_mod - av_mod, color='gray', alpha=0.5)
ax1.plot(inv_w, ((sac + ssac) / st[kvis] - 1.0) * av_mod, color='brown', lw=2, label='aC')
ax1.plot(inv_w, ((sasi + sssi) / st[kvis] - 1.0) * av_mod, color='orange', lw=2, label='Si')
ax1.plot(inv_w, extinction, color='magenta', lw=3, label='total')
ax1.set_xlim(-0.01, 11)
ax1.set_ylim(-2, 10.5)
ax1.set_xlabel(r'1/$\lambda$ ($\mu m^{-1}$)')
ax1.set_ylabel('E($\lambda$-V)')

# 2. Polarization Plot
wpol, spc, spsi, spd, stm = read_data('./Output/PolKappa.out')
wpol *= 1e4
spcd, spsid, spdd = [x * 100.0 / stm[kvis] for x in [spc, spsi, spd]]
sp = spdd + spcd + spsid

ax2 = axs[0, 1]
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.fill_between(wpol, 0.0021, spdd, color='lightgray')
ax2.plot(wpol, sp, color='magenta', lw=2)
ax2.plot(wpol, spsid, color='orange', lw=1.5)
ax2.set_xlim(0.3, 14)
ax2.set_ylim(0.0021, np.max(sp) * 1.1)
ax2.set_xlabel(r'$\lambda$ ($\mu m$)')
ax2.set_ylabel('p (%)')

# 3. Emission Plot
w_emi, emi_tot, emi_c, emi_si, _, emi_d, emi_vgr, emi_vsi, emipah, _ = read_data('./Output/emis.out', skiprows=2)
fak = 1.0 / np.max(emi_tot)

ax3 = axs[1, 0]
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.fill_between(w_emi, 4e-4, emi_d * fak, color='lightgray')
ax3.plot(w_emi, emi_tot * fak, color='magenta', lw=2)
ax3.set_xlim(2, 3000)
ax3.set_ylim(4e-4, 1.2)
ax3.set_xlabel(r'$\lambda$ ($\mu m$)')
ax3.set_ylabel(r'Intensity $I/I_{max}$')

# 4. Polarized Emission Plot
w_ep, emi_tot_p, emip_t, emip_c, emip_si, emip_d = read_data('./Output/emipol.out', skiprows=2)
fak_p = 1.0 / np.max(emi_tot_p)

ax4 = axs[1, 1]
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.plot(w_ep, emip_t * fak_p, color='blue', lw=2, label='total')
ax4.plot(w_ep, emip_si * fak_p, color='orange', lw=2, label='Si')
ax4.set_xlim(50, 3000)
ax4.set_ylim(1.01e-4, 1)
ax4.set_xlabel(r'$\lambda$ ($\mu m$)')
ax4.set_ylabel(r'Polarised Emission $P/I_{max}$')
ax4.legend()

plt.tight_layout()
plt.savefig('Example_RedPolEmis.pdf')
plt.show()
