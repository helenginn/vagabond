//
//  OpenGLView.m
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#import "OpenGLView.h"
#import "GLKeeper.h"

@interface OpenGLView ()
{
    GLKeeper *_glKeeper;
}

@property GLKeeper *glKeeper;
@property BOOL isSetup;
@property (nonatomic, retain) NSTimer *timer;
@property CGPoint lastMouse;

@end

@implementation OpenGLView

@synthesize glKeeper = _glKeeper;
@synthesize isSetup = _isSetup;
@synthesize timer = _timer;

- (void)awakeFromNib
{
    [self setupDepthBuffer];
}

- (void)setupGLKeeper {
    self.glKeeper = new GLKeeper(self.frame.size.width, self.frame.size.height);
    
    _isSetup = YES;
}

- (void)setupTimer
{
    _timer = [NSTimer scheduledTimerWithTimeInterval:1 / 20 target:self selector:@selector(render) userInfo:Nil repeats:YES];
}

- (void)prepareOpenGL {
    
    // this sets swap interval for double buffering
    GLint swapInt = 1;
    [[self openGLContext] setValues:&swapInt forParameter:NSOpenGLCPSwapInterval];
    
    // this enables alpha in the frame buffer (commented now)
    //  GLint opaque = 0;
    //  [[self openGLContext] setValues:&opaque forParameter:NSOpenGLCPSurfaceOpacity];
}

- (void)setWorkingDirectory
{
    NSString *resourcePath = [[NSBundle mainBundle] resourcePath];
    [[NSFileManager defaultManager] changeCurrentDirectoryPath:resourcePath];
}

- (void)setupFrameChange
{
    [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(changeFrame) name:@"NSViewFrameDidChangeNotification" object:nil];
}

- (void)mouseMoved:(NSEvent *)theEvent
{
    CGPoint newMouse = [theEvent locationInWindow];
    CGPoint movement = CGPointMake(newMouse.x - _lastMouse.x, newMouse.y - _lastMouse.y);
    
    self.glKeeper->movedMouse(movement.x, movement.y);
    
    _lastMouse = newMouse;
    
    if ((newMouse.x > self.frame.size.width - 20 || newMouse.x < 20) || (newMouse.y > self.frame.size.height - 20 || newMouse.y < 20))
    {
        CGPoint centre = CGPointMake(self.frame.size.width / 2, self.frame.size.height / 2);
        CGDisplayMoveCursorToPoint(kCGDirectMainDisplay, centre);
        _lastMouse = centre;
    }
    
    if (![_timer isValid])
        [self setupTimer];
}

- (void)mouseUp:(NSEvent *)theEvent
{
    [_timer invalidate];
}

-(void)keyDown:(NSEvent *)theEvent
{
    if ([theEvent modifierFlags] & NSNumericPadKeyMask)
    {
        [self interpretKeyEvents:[NSArray arrayWithObject:theEvent]];
    }
    
    NSString *characters = [theEvent charactersIgnoringModifiers];
    
    char character = [characters characterAtIndex:0];
    
    self.glKeeper->keyPressed(character);
    
    if (![_timer isValid])
        [self setupTimer];
        
}

- (void)keyUp:(NSEvent *)theEvent
{
    [_timer invalidate];
}

- (void)moveLeft:(id)sender
{
    self.glKeeper->rotateAngles(0.0, -0.02, 0);
}

- (void)moveUp:(id)sender
{
    self.glKeeper->rotateAngles(-0.02, 0.0, 0);
}

- (void)moveDown:(id)sender
{
    self.glKeeper->rotateAngles(0.02, 0.0, 0);
}

- (void)moveRight:(id)sender
{
    self.glKeeper->rotateAngles(0.0, 0.02, 0);
}

- (void)setup
{
    [self setupFrameChange];
    [self setWorkingDirectory];
    [self prepareOpenGL];
    [self setupGLKeeper];
    [self setupTimer];
    
    [self becomeFirstResponder];
    
    CGPoint centre = CGPointMake(self.frame.size.width / 2, self.frame.size.height / 2);
    CGDisplayMoveCursorToPoint(kCGDirectMainDisplay, centre);
    _lastMouse = centre;

    [NSCursor hide];
}

- (void)render
{
    [self setNeedsDisplay:YES];
}

- (void)setupDepthBuffer
{
    NSOpenGLPixelFormatAttribute attrs[] = {
        // NSOpenGLPFADoubleBuffer,
        NSOpenGLPFADepthSize, 32,
        0
    };
    NSOpenGLPixelFormat *format = [[NSOpenGLPixelFormat alloc] initWithAttributes:attrs];
    [self setPixelFormat:format];
}

- (void)drawRect:(NSRect)bounds {
    if (!_isSetup)
    {
        [self setup];
    }
    
    self.glKeeper->render();
    
    glFlush();
    
    int err;
    if ((err = glGetError()) != 0)
        NSLog(@"glGetError(): %d", err);
}

- (void)changeFrame
{
    self.glKeeper->changeSize(self.frame.size.width, self.frame.size.height);
}



@end
