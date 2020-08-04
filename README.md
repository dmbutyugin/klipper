This is a fork of [Klipper](https://www.klipper3d.org/)
([GitHub](https://github.com/KevinOConnor/klipper)), a 3d-Printer firmware.
The code in this repo adds experimental features to reduce resonances and
ringing during printing. The current focus is to extensively test adaptive
S-Curve acceleration - the input shaping feature was already integrated into
the mainline Klipper codebase.

The current branch recommended for the broad use is **`scurve-shaping`**
([installation and tuning instructions](https://github.com/dmbutyugin/klipper/blob/scurve-shaping/docs/Resonance_Compensation.md#switch-to-s-curve-acceleration-branch)).
All other branches are either deprecated in favor of this one, or have been
published for convenience and contain more experimental features and lack
proper documentation. Please update your installation to the suggested branch
if you are using one of the older branches.

The feedback is welcome at the
[Adaptive acceleration support](https://github.com/KevinOConnor/klipper/issues/3026)
feature request in the main Klipper repo. If you have issues specifically with
input shaping, please use
[this ticket](https://github.com/KevinOConnor/klipper/issues/3025)
instead. In general, if you report bugs, please double-check that the issue
cannot be reproduced on the mainline Klipper code without input shaping, and
attach the **full** `klippy.log` of the failed print attempt to your bug report.
Typically the problematic GCode is also required to debug the issue. However, if
the issue can be reproduced on the mainline Klipper code, please create a
separate issue
[here](https://github.com/KevinOConnor/klipper/issues) instead, and attach
`klippy.log` from the attempt on the **mainline** code.

Updates:

  * 2020-08-04: `[input_shaper]` is fully merged into the mainline Klipper
    (together with tuning docs).
  * 2020-07-09: `[input_shaper]` is merged into the mainline Klipper (except the
    docs for now)! Other approaches (e.g. `[smooth_axis]`) are incompatible with
    input shaping, which means that all S-Curve branches except `scurve-shaping`
    can no longer be updated. All users of `scurve-smoothing` and
   `scurve-c-combine-smoothpa` are encouraged to migrate to `scurve-shaping`
    branch (or to the mainline Klipper) and update their configurations
    [accordingly](https://github.com/dmbutyugin/klipper/blob/scurve-shaping/docs/Resonance_Compensation.md#switch-to-s-curve-acceleration-branch).
  * 2020-07-06: `scurve-shaping` branch is now recommended for broad use.
