// Real-time 2D Interference and Diffraction Simulator.
// by Umberto Puddu & Risa Charvi Metta
//
// Build:  make
// Run:    ./SlitDiffraction

#import <Cocoa/Cocoa.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

// Physics Configuration
static const double domainZ = 20.0;
static const double domainY = 10.0;
static const double ptsPerWave = 4.0;
static const int maxSources = 64;
static const double baseBrightness = 1.5;
static const double bgGlow = 0.1;
static const double fieldBlend = 0.8;
static const double glowBlend = 0.2;

// Resolution
static const int renderW = 480;
static const int renderH = 240;

// Slider Configurations
static const double lamMin = 380.0;
static const double lamMax = 780.0;
static const double lamInit = 500.0;

static const double wMin = 0.1;
static const double wMax = 5.0;
static const double wInit = 0.5;

static const double dMin = 1.0;
static const double dMax = 10.0;
static const double dInit = 3.0;

static const double nMin = 1.0;
static const double nMax = 5.0;
static const double nInit = 2.0;

static const double speedMin = 0.0;
static const double speedMax = 2.0;
static const double speedInit = 0.3;

// Physics Simulation
static void wavelengthToRGB(double lam_nm, double &r, double &g, double &b) {
  r = g = b = 0.0;

  // Clamp to visible light
  if (lam_nm < 440.0) {
    r = -(lam_nm - 440.0) / 60.0;
    b = 1.0;
  } else if (lam_nm < 490.0) {
    g = (lam_nm - 440.0) / 50.0;
    b = 1.0;
  } else if (lam_nm < 510.0) {
    g = 1.0;
    b = -(lam_nm - 510.0) / 20.0;
  } else if (lam_nm < 580.0) {
    r = (lam_nm - 510.0) / 70.0;
    g = 1.0;
  } else if (lam_nm < 645.0) {
    r = 1.0;
    g = -(lam_nm - 645.0) / 65.0;
  } else {
    r = 1.0;
  }

  // Clamp invalid colors
  if (r < 0)
    r = 0;
  if (r > 1)
    r = 1;
  if (g < 0)
    g = 0;
  if (g > 1)
    g = 1;
  if (b < 0)
    b = 0;
  if (b > 1)
    b = 1;
}

static std::vector<double> simulateWaves(double lam_nm, double slit_w_um,
                                         double slit_sep_um, int n_slits,
                                         int nx, int ny, double Z_um,
                                         double Y_um, double time_phase) {
  const double lam = lam_nm * 1e-3; // nm to um
  const double w = slit_w_um;
  const double d = slit_sep_um;
  const double k = 2.0 * M_PI / lam; // Wave number

  int src_per_slit = std::max(1, (int)(w / (lam / ptsPerWave)));
  if (src_per_slit > maxSources)
    src_per_slit = maxSources;

  std::vector<double> sources;
  sources.reserve(n_slits * src_per_slit);
  for (int s = 0; s < n_slits; ++s) {
    const double centre = (n_slits == 1) ? 0.0
                                         : (s - 0.5 * (n_slits - 1)) *
                                               d; // Find center of this slit
    for (int i = 0; i < src_per_slit; ++i) {
      const double t = (i + 0.5) / src_per_slit -
                       0.5; // Distribute points evenly across the slit width
      sources.push_back(centre + w * t);
    }
  }

  std::vector<double> field(nx * ny);
  double *field_ptr = field.data();

  dispatch_apply(
      ny, dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0),
      ^(size_t py) {
        const double y = -Y_um * 0.5 + (py + 0.5) * (Y_um / ny);
        for (int px = 0; px < nx; ++px) {
          const double z = (px + 0.5) * (Z_um / nx);
          double E = 0.0;
          // For z very close to 0, limit 1/sqrt(r) (to avoid NaN)
          double min_r = lam / ptsPerWave;
          for (double ys : sources) {
            const double dy = y - ys;
            double r = std::sqrt(
                z * z + dy * dy); // Distance from point source to pixel
            if (r < min_r)
              r = min_r;
            const double phi = k * r - time_phase; // Instantaneous phase
            E += std::cos(phi) /
                 std::sqrt(r); // Superposition of 2D cylindrical wave
          }
          field_ptr[py * nx + px] = E;
        }
      });
  return field;
}

static NSImage *renderFrame(double lam_nm, double slit_w_um, double slit_sep_um,
                            int n_slits, int width_px, int height_px,
                            double time_phase) {
  auto field = simulateWaves(lam_nm, slit_w_um, slit_sep_um, n_slits, width_px,
                             height_px, domainZ, domainY, time_phase);

  double cr, cg, cb;
  wavelengthToRGB(lam_nm, cr, cg, cb);

  // Estimate of maximum field amplitude
  double expected_max =
      (n_slits * std::max(1.0, slit_w_um / (lam_nm * 1e-3 / ptsPerWave))) /
      std::sqrt(domainZ * 0.1);
  if (expected_max <= 0)
    expected_max = 1.0;

  double scale = baseBrightness / expected_max;

  NSBitmapImageRep *rep =
      [[NSBitmapImageRep alloc] initWithBitmapDataPlanes:NULL
                                              pixelsWide:width_px
                                              pixelsHigh:height_px
                                           bitsPerSample:8
                                         samplesPerPixel:3
                                                hasAlpha:NO
                                                isPlanar:NO
                                          colorSpaceName:NSDeviceRGBColorSpace
                                             bytesPerRow:width_px * 3
                                            bitsPerPixel:24];

  unsigned char *buf = [rep bitmapData];
  double *field_ptr = field.data();

  dispatch_apply(height_px,
                 dispatch_get_global_queue(QOS_CLASS_USER_INTERACTIVE, 0),
                 ^(size_t py) {
                   unsigned char *row = buf + py * (width_px * 3);
                   for (int px = 0; px < width_px; ++px) {
                     double E = field_ptr[py * width_px + px];
                     double val = E * scale;

                     // Draw only positive crests
                     double v = val > 0 ? val : 0.0;
                     if (v > 1.0)
                       v = 1.0;

                     // Glow to highlight fringes
                     double intensity = E * E * scale * scale * bgGlow;
                     if (intensity > 1.0)
                       intensity = 1.0;

                     v = v * fieldBlend + intensity * glowBlend;
                     if (v > 1.0)
                       v = 1.0;

                     row[3 * px + 0] = (unsigned char)(v * cr * 255.0);
                     row[3 * px + 1] = (unsigned char)(v * cg * 255.0);
                     row[3 * px + 2] = (unsigned char)(v * cb * 255.0);
                   }
                 });

  NSImage *img = [[NSImage alloc] initWithSize:NSMakeSize(width_px, height_px)];
  [img addRepresentation:rep];
  return img;
}

@interface AppDelegate : NSObject <NSApplicationDelegate>
@property(strong) NSWindow *window;
@property(strong) NSImageView *imageView;
@property(strong) NSSlider *lamSlider;
@property(strong) NSSlider *wSlider;
@property(strong) NSSlider *dSlider;
@property(strong) NSSlider *nSlider;
@property(strong) NSSlider *speedSlider;
@property(strong) NSTextField *lamLabel;
@property(strong) NSTextField *wLabel;
@property(strong) NSTextField *dLabel;
@property(strong) NSTextField *nLabel;
@property(strong) NSTextField *speedLabel;
@property(strong) NSTextField *statusLabel;
@property(strong) NSButton *playButton;
@property(strong) NSTimer *timer;
@property(assign) double time_phase;
@end

@implementation AppDelegate

- (NSSlider *)createSlider:(CGFloat)y
                  minValue:(double)mn
                  maxValue:(double)mx
                   initial:(double)init
                      step:(double)step {
  NSSlider *s = [[NSSlider alloc] initWithFrame:NSMakeRect(140, y, 740, 22)];
  s.minValue = mn;
  s.maxValue = mx;
  s.doubleValue = init;
  if (step > 0) {
    s.allowsTickMarkValuesOnly = YES;
    s.numberOfTickMarks = (int)((mx - mn) / step) + 1;
  }
  s.target = self;
  s.action = @selector(onSliderChanged:);
  s.continuous = YES;
  [self.window.contentView addSubview:s];
  return s;
}

- (NSTextField *)createLabel:(CGFloat)y text:(NSString *)t width:(CGFloat)w {
  NSTextField *f =
      [[NSTextField alloc] initWithFrame:NSMakeRect(20, y - 2, w, 22)];
  f.stringValue = t;
  f.bezeled = NO;
  f.editable = NO;
  f.selectable = NO;
  f.drawsBackground = NO;
  f.font = [NSFont systemFontOfSize:13 weight:NSFontWeightMedium];
  [self.window.contentView addSubview:f];
  return f;
}

- (NSTextField *)createReadout:(CGFloat)y {
  NSTextField *f =
      [[NSTextField alloc] initWithFrame:NSMakeRect(890, y - 2, 90, 22)];
  f.bezeled = NO;
  f.editable = NO;
  f.selectable = NO;
  f.drawsBackground = NO;
  f.alignment = NSTextAlignmentRight;
  f.font = [NSFont monospacedDigitSystemFontOfSize:13
                                            weight:NSFontWeightRegular];
  [self.window.contentView addSubview:f];
  return f;
}

- (void)onLaunch {
  NSRect frame = NSMakeRect(120, 120, 1000, 750);
  self.window = [[NSWindow alloc]
      initWithContentRect:frame
                styleMask:(NSWindowStyleMaskTitled | NSWindowStyleMaskClosable |
                           NSWindowStyleMaskMiniaturizable |
                           NSWindowStyleMaskFullSizeContentView)
                  backing:NSBackingStoreBuffered
                    defer:NO];
  self.window.titlebarAppearsTransparent = YES;
  self.window.title = @"Slit Diffraction";
  [self.window center];

  NSVisualEffectView *visualEffectView =
      [[NSVisualEffectView alloc] initWithFrame:frame];
  visualEffectView.material = NSVisualEffectMaterialUnderWindowBackground;
  visualEffectView.blendingMode = NSVisualEffectBlendingModeBehindWindow;
  visualEffectView.state = NSVisualEffectStateActive;
  self.window.contentView = visualEffectView;

  self.imageView =
      [[NSImageView alloc] initWithFrame:NSMakeRect(20, 230, 960, 480)];
  self.imageView.imageScaling = NSImageScaleAxesIndependently;
  self.imageView.wantsLayer = YES;
  self.imageView.layer.cornerRadius = 8.0;
  self.imageView.layer.masksToBounds = YES;
  [self.window.contentView addSubview:self.imageView];

  [self createLabel:205 text:@"λ (nm)" width:120];
  self.lamSlider = [self createSlider:205
                             minValue:lamMin
                             maxValue:lamMax
                              initial:lamInit
                                 step:0];
  self.lamLabel = [self createReadout:205];

  [self createLabel:170 text:@"Width (μm)" width:120];
  self.wSlider = [self createSlider:170
                           minValue:wMin
                           maxValue:wMax
                            initial:wInit
                               step:0];
  self.wLabel = [self createReadout:170];

  [self createLabel:135 text:@"Slit Distance (μm)" width:130];
  self.dSlider = [self createSlider:135
                           minValue:dMin
                           maxValue:dMax
                            initial:dInit
                               step:0];
  self.dLabel = [self createReadout:135];

  [self createLabel:100 text:@"Slits" width:120];
  self.nSlider = [self createSlider:100
                           minValue:nMin
                           maxValue:nMax
                            initial:nInit
                               step:1];
  self.nLabel = [self createReadout:100];

  [self createLabel:65 text:@"Wave Speed" width:120];
  self.speedSlider = [self createSlider:65
                               minValue:speedMin
                               maxValue:speedMax
                                initial:speedInit
                                   step:0];
  self.speedLabel = [self createReadout:65];

  self.playButton = [[NSButton alloc] initWithFrame:NSMakeRect(20, 20, 80, 28)];
  [self.playButton setButtonType:NSButtonTypePushOnPushOff];
  self.playButton.title = @"Play";
  self.playButton.bezelStyle = NSBezelStyleRounded;
  self.playButton.target = self;
  self.playButton.action = @selector(onTogglePlay:);
  self.playButton.font = [NSFont systemFontOfSize:13 weight:NSFontWeightMedium];
  [self.window.contentView addSubview:self.playButton];

  self.statusLabel =
      [[NSTextField alloc] initWithFrame:NSMakeRect(120, 24, 860, 22)];
  self.statusLabel.bezeled = NO;
  self.statusLabel.editable = NO;
  self.statusLabel.selectable = NO;
  self.statusLabel.drawsBackground = NO;
  self.statusLabel.textColor = [NSColor secondaryLabelColor];
  self.statusLabel.font = [NSFont systemFontOfSize:12
                                            weight:NSFontWeightRegular];
  [self.window.contentView addSubview:self.statusLabel];

  self.time_phase = 0.0;

  [self.window makeKeyAndOrderFront:nil];
  [NSApp activateIgnoringOtherApps:YES];
  [self onSliderChanged:nil];

  [[NSNotificationCenter defaultCenter]
      addObserver:self
         selector:@selector(onClose:)
             name:NSWindowWillCloseNotification
           object:self.window];
}

- (void)onClose:(NSNotification *)note {
  [NSApp terminate:nil];
}

- (void)onTogglePlay:(id)sender {
  if (self.playButton.state == NSControlStateValueOn) {
    self.playButton.title = @"Pause";
    self.timer = [NSTimer scheduledTimerWithTimeInterval:1.0 / 60.0
                                                  target:self
                                                selector:@selector(onTimerTick:)
                                                userInfo:nil
                                                 repeats:YES];
    [[NSRunLoop currentRunLoop] addTimer:self.timer
                                 forMode:NSRunLoopCommonModes];
  } else {
    self.playButton.title = @"Play";
    [self.timer invalidate];
    self.timer = nil;
  }
}

- (void)onTimerTick:(NSTimer *)t {
  self.time_phase += self.speedSlider.doubleValue;
  [self onSliderChanged:nil];
}

- (void)onSliderChanged:(id)sender {
  double lam = self.lamSlider.doubleValue;
  double w = self.wSlider.doubleValue;
  double d = self.dSlider.doubleValue;
  int n = (int)self.nSlider.intValue;

  self.imageView.image =
      renderFrame(lam, w, d, n, renderW, renderH, self.time_phase);
  self.lamLabel.stringValue = [NSString stringWithFormat:@"%6.1f", lam];
  self.wLabel.stringValue = [NSString stringWithFormat:@"%6.2f", w];
  self.dLabel.stringValue = [NSString stringWithFormat:@"%6.2f", d];
  self.nLabel.stringValue = [NSString stringWithFormat:@"%6d", n];
  self.speedLabel.stringValue =
      [NSString stringWithFormat:@"%6.2fx", self.speedSlider.doubleValue];

  if (n >= 2) {
    self.statusLabel.stringValue = [NSString
        stringWithFormat:@"Z = %.1f μm, Y = %.1f μm. %d-slit interference.",
                         domainZ, domainY, n];
  } else {
    self.statusLabel.stringValue = [NSString
        stringWithFormat:@"Z = %.1f μm, Y = %.1f μm. 1-slit diffraction.",
                         domainZ, domainY];
  }
}

@end

int main(int, const char *[]) {
  @autoreleasepool {
    [NSApplication sharedApplication];
    AppDelegate *del = [[AppDelegate alloc] init];
    [NSApp setDelegate:del];
    [NSApp setActivationPolicy:NSApplicationActivationPolicyRegular];
    [del onLaunch];
    [NSApp run];
    (void)del;
  }
  return 0;
}
