//
//  OpenGLView.m
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#import "OpenGLView.h"
#import "GLKeeper.h"

#define SLOW_FRAME_RATE 0.5

@interface OpenGLView ()
{
	GLKeeper *_glKeeper;
}

@property GLKeeper *glKeeper;
@property BOOL isSetup;
@property (nonatomic, retain) NSTimer *timer;
@property CGPoint lastMouse;
@property BOOL optionPressed;

@end

@implementation OpenGLView

@synthesize glKeeper = _glKeeper;
@synthesize isSetup = _isSetup;
@synthesize timer = _timer;
@synthesize optionPressed = _optionPressed;

- (void)awakeFromNib
{
	[self setupDepthBuffer];
}

- (void)setupGLKeeper {
	self.glKeeper = new GLKeeper(self.frame.size.width, self.frame.size.height);

	_isSetup = YES;
}

- (void)setupTimer:(double)interval
{
	_timer = [NSTimer scheduledTimerWithTimeInterval:interval target:self selector:@selector(render) userInfo:Nil repeats:YES];
}

- (void)prepareOpenGL {

	// this sets swap interval for double buffering
	GLint swapInt = 1;
	[[self openGLContext] setValues:&swapInt forParameter:NSOpenGLCPSwapInterval];

	CGFloat scale = [[NSScreen mainScreen] backingScaleFactor];

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

- (void)rightMouseDragged:(NSEvent *)event
{
	[self aMouseDragged:event];
}

- (void)mouseDragged:(NSEvent *)theEvent
{
	[self aMouseDragged:theEvent];
}

- (void)aMouseDragged:(NSEvent *)theEvent
{
	[self setupTimer:0.2];
	CGPoint newMouse = [theEvent locationInWindow];
	CGPoint movement = CGPointMake(newMouse.x - _lastMouse.x, newMouse.y - _lastMouse.y);

	if ([theEvent type] == NSEventTypeLeftMouseDragged && !_optionPressed)
	{
		self.glKeeper->draggedLeftMouse(movement.x, movement.y);
	}
	else if ([theEvent type] == NSEventTypeLeftMouseDragged && _optionPressed)
	{
		self.glKeeper->panned(movement.x, movement.y);
	}
	else if ([theEvent type] == NSEventTypeRightMouseDragged)
	{
		self.glKeeper->draggedRightMouse(movement.x, movement.y);
	}

	_lastMouse = newMouse;

	if (![_timer isValid])
		[self setupTimer:(0.2)];
}

- (void)rightMouseDown:(NSEvent *)theEvent
{
	[self aMouseDown:theEvent];
}

- (void)mouseDown:(NSEvent *)theEvent
{
	[self aMouseDown:theEvent];
}

- (void)aMouseDown:(NSEvent *)theEvent
{
	_lastMouse = [theEvent locationInWindow];
	[self setupTimer:SLOW_FRAME_RATE];
}

-(void)flagsChanged:(NSEvent *)theEvent
{
	if ([theEvent modifierFlags] & NSEventModifierFlagControl)
	{
		_optionPressed = TRUE;
	}
	else
	{
		_optionPressed = FALSE;
	}
}

- (void)keyUp:(NSEvent *)theEvent
{
	_optionPressed = FALSE;
	[self setupTimer:(SLOW_FRAME_RATE)];

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
	[self setWorkingDirectory];
	[self prepareOpenGL];
	[self setupGLKeeper];
	[self setupTimer:(0.2)];
	[self setupFrameChange];

	[self becomeFirstResponder];

	CGPoint centre = CGPointMake(self.frame.size.width / 2, self.frame.size.height / 2);
	_lastMouse = centre;
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

- (void)drawRect:(NSRect)bounds
{
	if (!_isSetup)
	{
		[self setup];
	}

	CGFloat scale = [[NSScreen mainScreen] backingScaleFactor];
	CALayer *eaglLayer = [self layer];
	eaglLayer.contentsScale = scale;

	CGSize s = self.frame.size;
	glViewport(0, 0, s.width * scale, s.height * scale);
	_glKeeper->changeSize(s.width, s.height);

	self.glKeeper->render();
	glFlush();

	int err;
	if ((err = glGetError()) != 0)
		NSLog(@"glGetError(): %d", err);
}

- (void)changeFrame
{
	return;

	CGFloat scale = [[NSScreen mainScreen] backingScaleFactor];
	glViewport(0, 0, self.frame.size.width * scale,
			   self.frame.size.height * scale);

}



@end
