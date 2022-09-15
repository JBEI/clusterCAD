import React from 'react';
import { connect } from 'react-redux';

function mapStateToProps(state) {
  const { domainSearchResponse } = state;
  return { responseData: domainSearchResponse };
}

class SearchResults extends React.Component {

  constructor(props) {
    super(props);
  }

  render() {
    return(
      <iframe className="Results form" srcdoc={this.props.responseData}>
        
      </iframe>
    )
  }

}

export default connect(mapStateToProps)(SearchResults);