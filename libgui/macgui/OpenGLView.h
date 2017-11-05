//
//  OpenGLView.h
//  RaddoseViewer
//
//  Created by Helen Ginn on 19/05/2014.
//  Copyright (c) 2014 Raddose3D. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <QuartzCore/QuartzCore.h>
#include <OpenGL/OpenGL.h>


@interface OpenGLView : NSOpenGLView
{
    
}

- (void) drawRect: (NSRect) bounds;
-(IBAction)keyDown:(NSEvent *)theEvent;
- (void)setup;

@end
