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
      <div className="SearchResults">
        <iframe className="Results form" srcDoc={this.props.responseData}>
          
        </iframe>
      </div>
    )
  }

}

export default connect(mapStateToProps)(SearchResults);