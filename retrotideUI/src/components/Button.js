import React from 'react';

class Button extends React.Component {
  render() {
    let {href, disabled, className, onClick, children} = this.props;
    // we need to put the className here so the styling is inclusive of the a tag
    // note that the Boolean 'disabled' value MUST be in curly braces and not quotes
    // when it is passed in
    return(
      <button href={href} disabled={disabled} className={className + ' button'} onClick={onClick} >
        {children}
      </button>
    )
  }
}

export default Button;